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
	public:
		virtual void SyncTo(Oscillator& dst)
		{
			syncDst = &dst;
		}
		virtual void UnregSync()
		{
			syncDst = nullptr;
		}

		virtual float GetDT() = 0;//获取目前的采样增量
		virtual bool IsWrapThisSample() = 0;//这个采样是否发生wrap
		virtual float GetWrapWhere() = 0;//获取目前的采样在哪里wrap

		virtual void Step(float dt) = 0;//分开是为了更好处理振荡器内复杂状态
		virtual float Get() = 0;

	};

	class SawOscillator :public Oscillator
	{
	private:
		TableBlep blep;
		float t = 0;
		float dt = 0;

		int isWrap = 0;
		float t2 = 0, amp2 = 0, where2 = 0;

		float startPhase = 0;//初始相位 [0,1]
	public:
		void SetStartPhase(float phi)
		{
			startPhase = phi;
		}
		void SetPhase(float phi, float where = 0)
		{
			float newt = phi;
			blep.Add(newt - t, where);
			t = newt;
		}

		bool IsWrapThisSample() override
		{
			return isWrap;
		}
		float GetWrapWhere() override
		{
			return where2;
		}
		float GetDT() override
		{
			return dt;
		}
		void Step(float dt1) override
		{
			dt = dt1;
			if (dt > 1.0)dt = 1.0;
			if (dt < -1.0)dt = -1.0;
			t += dt;

			isWrap = 0;
			if (t >= 1.0)WrapPhaseDown();
			else if (t < 0.0)WrapPhaseUp();
		}

		void WrapPhaseDown()//只预备状态，不更新
		{
			t2 = t;
			amp2 = -1.0;
			t2 -= 1.0;
			where2 = t2 / dt;
			if (where2 < 0.0)where2 = 0.0;
			if (where2 > 1.0)where2 = 1.0;
			isWrap = 1;
		}
		void WrapPhaseUp()//只预备状态，不更新
		{
			t2 = t;
			where2 = t2 / dt;
			if (where2 < 0.0)where2 = 0.0;
			if (where2 > 1.0)where2 = 1.0;
			amp2 = 1.0;
			t2 += 1.0;
			isWrap = 1;
		}
		void ApplyWrap()
		{
			t = t2;
			blep.Add(amp2, where2);
		}
		void DoSync(float dstWhere)
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
		float Get() override
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

	class WaveformOsc :public Oscillator
	{
	private:
		SawOscillator saw1, saw2;
		float duty = 0.25, dt = 0.0;
		float tri = 0;
		float form = 0;//imp(0)->pwm(1)->tri(2)

		float lastv = 0;
		float lastpwm = 0;
	public:
		WaveformOsc()
		{
			SetPWM(0.25);
		}
		void SetPWM(float duty)
		{
			duty *= 0.5;
			this->duty = duty;
			saw1.SetStartPhase(0);
			saw2.SetStartPhase(duty);
		}
		void SetWaveform(float form)
		{
			this->form = form;
		}
		void SyncTo(Oscillator& dst) override
		{
			Oscillator::SyncTo(dst);
			saw1.SyncTo(dst);
			saw2.SyncTo(dst);
		}
		void UnregSync() override
		{
			Oscillator::UnregSync();
			saw1.UnregSync();
			saw2.UnregSync();
		}

		bool IsWrapThisSample() override
		{
			return saw1.IsWrapThisSample();
		}
		float GetWrapWhere() override
		{
			return saw1.GetWrapWhere();
		}
		float GetDT() override
		{
			return dt;
		}
		void Step(float _dt) override
		{
			dt = _dt;
			if (dt > 1.0)dt = 1.0;
			if (dt < -1.0)dt = -1.0;
			saw1.Step(dt);
			saw2.Step(dt);
		}
		float Get() override
		{
			float v1 = saw1.Get();
			float v2 = saw2.Get();
			float pwm = v1 - v2;

			float pwmdc = 0;
			if (syncDst)
			{
				float syncdt = syncDst->GetDT();
				float k = dt / syncdt;
				float r = k - (int)k;
				if (r <= 1.0f - duty) {
					pwmdc = -2.0f * duty * r / k;
				}
				else {
					pwmdc = 2.0f * (1.0f - duty) * (r - 1.0f) / k;
				}
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

	class OscTest
	{
	private:
		WaveformOsc osc1;
		WaveformOsc osc2;
		float dt1 = 0, dt2 = 0;
		float duty = 0.5;
	public:
		void SetParams(float freq, float sync, float pwm, float form, float sr)
		{
			dt1 = freq / sr;
			dt2 = freq * sync / sr;
			this->duty = pwm;
			osc2.SyncTo(osc1);
			osc2.SetPWM(duty);
			osc2.SetWaveform(form);
		}
		float ProcessSample()
		{
			osc1.Step(dt1);
			osc2.Step(dt2);
			float v1 = osc1.Get();
			float v2 = osc2.Get();
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
	};
}