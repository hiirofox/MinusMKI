#pragma once

#include "TableBlep.h"
#include <functional>
#include <cassert>

namespace MinusMKI
{
	class Oscillator
	{
	private:
	protected:
		Oscillator* syncDst = nullptr;
		using Blep = TableBlep;
	public:
		virtual void SyncTo(Oscillator& dst)
		{
			syncDst = &dst;
		}
		virtual void UnregSync()
		{
			syncDst = nullptr;
		}

		virtual float GetDT() const = 0;//获取目前的采样增量
		virtual bool IsWrapThisSample() const = 0;//这个采样是否发生wrap
		virtual float GetWrapWhere() const = 0;//获取目前的采样在哪里wrap

		virtual void Step(float dt) = 0;//分开是为了更好处理振荡器内复杂状态
		virtual float Get() = 0;

	};

	class SawOscillator final :public Oscillator
	{
	private:
		Blep blep;
		float t = 0;
		float dt = 0;

		int isWrap = 0;
		float t2 = 0, amp2 = 0, where2 = 0;

		float startPhase = 0;//初始相位 [0,1]
	public:
		inline void SetStartPhase(float phi)
		{
			startPhase = phi;
		}
		inline void SetPhase(float phi, float where = 0)
		{
			float newt = phi;
			blep.Add(newt - t, where);
			t = newt;
		}

		inline bool IsWrapThisSample() const final override
		{
			return isWrap;
		}
		inline float GetWrapWhere() const  final override
		{
			return where2;
		}
		inline float GetDT() const final override
		{
			return dt;
		}
		inline void Step(float dt1) final override
		{
			dt = dt1;
			if (dt > 1.0)dt = 1.0;
			if (dt < -1.0)dt = -1.0;
			t += dt;

			isWrap = 0;
			if (t >= 1.0)WrapPhaseDown();
			else if (t < 0.0)WrapPhaseUp();
		}

		inline void WrapPhaseDown()//只预备状态，不更新
		{
			t2 = t;
			amp2 = -1.0;
			t2 -= 1.0;
			where2 = t2 / dt;
			if (where2 < 0.0)where2 = 0.0;
			if (where2 > 1.0)where2 = 1.0;
			isWrap = 1;
		}
		inline void WrapPhaseUp()//只预备状态，不更新
		{
			t2 = t;
			where2 = t2 / dt;
			if (where2 < 0.0)where2 = 0.0;
			if (where2 > 1.0)where2 = 1.0;
			amp2 = 1.0;
			t2 += 1.0;
			isWrap = 1;
		}
		inline void ApplyWrap()
		{
			t = t2;
			blep.Add(amp2, where2);
		}
		inline void DoSync(float dstWhere)
		{
			if (isWrap && dstWhere < where2)
			{
				ApplyWrap();//先处理自己的wrap
				isWrap = 0;
			}
			float syncPhase = dstWhere * dt + startPhase;
			float diff = syncPhase - t;
			blep.Add(diff, dstWhere);
			t = syncPhase;
			if (t >= 1.0)
			{
				int amp = t;
				t -= amp;
				blep.Add(-amp, t / dt);
			}
		}
		float Get() final override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncTime = isSync ? syncDst->GetWrapWhere() : -1.0f;

			if (isSync)
			{
				DoSync(syncTime);
			}
			else if (isWrap)
			{
				ApplyWrap();
			}

			blep.Step();
			float v = t + blep.Get();
			return v * 2.0 - 1.0;
		}
	};

	class WaveformOsc final :public Oscillator
	{
	private:
		SawOscillator saw1, saw2;
		float duty = 0.25, dt = 0.0;
		float tri = 0;
		float form = 0;//imp(0)->pwm(1)->tri(2)

		float lastv = 0;
		float lastpwm = 0;
		float startPhase = 0;
	public:
		WaveformOsc()
		{
			UnregSync();
			SetPWM(0.25);
		}
		inline void SetStartPhase(float sp)
		{
			startPhase = sp;
			float p1 = startPhase;
			float p2 = startPhase + duty;
			p1 -= (int)p1;
			p2 -= (int)p2;
			saw1.SetStartPhase(p1);
			saw2.SetStartPhase(p2);
		}
		inline void SetPWM(float duty)
		{
			duty *= 0.5;
			this->duty = duty;
			float p1 = startPhase;
			float p2 = startPhase + duty;
			p1 -= (int)p1;
			p2 -= (int)p2;
			saw1.SetStartPhase(p1);
			saw2.SetStartPhase(p2);
		}
		inline void SetWaveform(float form)
		{
			this->form = form;
		}
		void SyncTo(Oscillator& dst) final override
		{
			Oscillator::SyncTo(dst);
			saw1.SyncTo(dst);
			saw2.SyncTo(dst);
		}
		void UnregSync() final override
		{
			/*
			Oscillator::UnregSync();
			saw1.UnregSync();
			saw2.UnregSync();
			*/
			SyncTo(*this);//你别管他，去掉就跑不了了
		}

		inline bool IsWrapThisSample() const final override
		{
			return saw1.IsWrapThisSample();
		}
		inline float GetWrapWhere() const final override
		{
			return saw1.GetWrapWhere();
		}
		inline float GetDT() const final override
		{
			return dt;
		}
		inline void Step(float _dt) final override
		{
			dt = _dt;
			if (dt > 1.0)dt = 1.0;
			if (dt < -1.0)dt = -1.0;
			saw1.Step(dt);
			saw2.Step(dt);
		}
		inline float GetPwmIntegral(float p, float d) const {
			p -= (int)p;
			if (p < 0.0f) p += 1.0f;
			if (p <= 1.0f - d) {
				return -2.0f * d * p;
			}
			else {
				return 2.0f * (1.0f - d) * (p - 1.0f);
			}
		}
		float Get() final override
		{
			float v1 = saw1.Get();
			float v2 = saw2.Get();
			float pwm = v1 - v2;
			//return pwm;//test
			
			float pwmdc = 0;
			if (syncDst)
			{
				float syncdt = syncDst->GetDT();
				float k = dt / syncdt;
				float phaseStart = startPhase;
				float phaseEnd = startPhase + k;
				float intStart = GetPwmIntegral(phaseStart, duty);
				float intEnd = GetPwmIntegral(phaseEnd, duty);
				pwmdc = (intEnd - intStart) / k;
			}
			pwm -= pwmdc; // 减去计算出的 DC，得到纯正的无交流偏移波形

			float dutyfix = duty + 0.001;
			float trifix = dt / (dutyfix * (1.0f - dutyfix));
			tri = pwm * trifix + tri * 0.999;

			float fixmix = duty * 50;
			if (fixmix <= 1.0)
			{
				tri = -v1 * (1.0f - fixmix) + tri * fixmix;
			}

			if (fixmix <= 1.0)
			{
				float imp2 = v1 - lastv;
				pwm = -imp2 * (1.0f - fixmix) + pwm * fixmix;
			}
			lastv = v1;

			float imp = pwm - lastpwm;
			lastpwm = pwm;

			float mix1 = 0, mix2 = 0, mix3 = 0;
			if (form <= 1.0f)
			{
				mix1 = 1.0f - form;
				mix2 = form;
			}
			else
			{
				mix2 = 2.0f - form;
				mix3 = form - 1.0f;
			}

			float out = imp * mix1 + pwm * mix2 + tri * mix3;
			return out;
		}
	};

	class ZeroCrossDetector final :public Oscillator
	{
	private:
		float z0 = 0;
		float z1 = 0;
		float z2 = 0;
		float z3 = 0;
		int isWrap = 0;
		float where = 0;

		float dt = 0.0000001;
		float cycleLen = 0;
	public:
		void SyncTo(Oscillator& dst) final override
		{
		}
		void UnregSync() final override
		{
		}

		inline bool IsWrapThisSample() const final override
		{
			if (syncDst) return syncDst->IsWrapThisSample();
			return isWrap;
		}
		inline float GetWrapWhere() const final override
		{
			if (syncDst) return syncDst->GetWrapWhere();
			return where;
		}
		inline float GetDT() const final override
		{
			if (syncDst) return syncDst->GetDT();
			return dt;
		}

		inline float hornor_fdivdf(float x, float c0, float c1, float c2, float c3)
		{
			float fx = c0 + x * (c1 + x * (c2 + x * c3));
			float dfx = c1 + x * (2.0f * c2 + x * 3.0f * c3);
			return fx / dfx;
		}
		inline void Step(float x) final override
		{
			z3 = z2;
			z2 = z1;
			z1 = z0;
			z0 = x;
			isWrap = 0;
			cycleLen += 1.0;
			//if (!(z2 <= 0.0f && z1 > 0.0f))return;//fake zero crossing
			//true zero crossing:

			float a = (z3 + z2) * 0.5;
			float b = (z1 + z0) * 0.5;
			//float a = z2;
			//float b = z1;
			float x0 = a / (a - b);

			//if (x0 < 0)x0 = 0;
			//if (x0 > 1)x0 = 1;
			if (x0 >= 0.0f && x0 <= 1.0f)
			{
				isWrap = 1;
				where = x0;
				dt = 1.0 / cycleLen;
				cycleLen = where;
			}
			//else assert(0);
		}
		float Get() final override
		{
			return 0;
		}
	};

	class OscTest
	{
	private:
		ZeroCrossDetector zcd;
		WaveformOsc osc1;
		WaveformOsc osc2;
		float dt1 = 0, dt2 = 0;
		float duty = 0.5;
		float fb = 0, fbv = 0;
	public:
		void SetParams(float freq, float sync, float pwm, float form, float fb, float sr)
		{
			dt1 = freq / sr;
			dt2 = freq * sync / sr;
			this->fb = fb;
			this->duty = pwm;
			osc1.SetPWM(duty);
			osc1.SetWaveform(form);

			osc2.SyncTo(osc1);
			osc2.SetPWM(duty);
			osc2.SetWaveform(form);

		}
		float lastv1 = 0;
		float ProcessSample()
		{
			zcd.Step(lastv1);
			osc1.Step(dt1 + fbv * fb);
			osc2.Step(dt2 + fbv * fb);
			float v1 = osc1.Get();
			float v2 = osc2.Get();
			lastv1 = v1;
			fbv = v2;
			return v2;
		}
		void ProcessBlock(float* outl, float* outr, int numSamples)
		{
			for (int i = 0; i < numSamples; ++i)
			{
				float v = ProcessSample() * 0.125;
				outl[i] = v;
				outr[i] = v;
			}
		}
		void SetStartPhase(float sp)
		{
			osc1.SetStartPhase(sp);
			osc2.SetStartPhase(sp);
		}
	};

	class UnisonTest2
	{
	private:
		constexpr static int UnisonNum = 1;
		OscTest wav[UnisonNum];
		float unitvol = 1.0 / sqrtf(UnisonNum);
	public:
		UnisonTest2()
		{
			for (int i = 0; i < UnisonNum; ++i)
			{
				float randphase = (float)rand() / RAND_MAX;
				wav[i].SetStartPhase(randphase);
			}
		}
		void SetParams(float freq, float sync, float pwm, float form, float fb, float detune, float sr)
		{
			for (int i = 0; i < UnisonNum; ++i)
			{
				float f = freq * (1.0 + ((float)i / UnisonNum - 0.5) * detune * 0.05);
				wav[i].SetParams(f, sync, pwm, form, fb, sr);
			}
		}
		void ProcessBlock(float* outl, float* outr, int numSamples)
		{
			for (int i = 0; i < numSamples; ++i)
			{
				float vl = ProcessSample();
				float vr = vl;

				outl[i] = vl / 8.0;
				outr[i] = vr / 8.0;
			}
		}
		inline float ProcessSample()
		{
			float sum = 0;
			for (int i = 0; i < UnisonNum; ++i)
			{
				sum += wav[i].ProcessSample();
			}
			sum *= unitvol;
			return sum;
		}
	};
}