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
	class TriOscillator final : public Oscillator
	{
	private:
		Blep blep;
		float t = 0;
		float dt = 0;

		// 周期边界 Wrap 状态 (底部峰值, 相位 0 或 1)
		int isWrap = 0;
		float t2 = 0, where2 = 0;

		// 内部占空比跨越状态 (顶部峰值, 相位 = duty)
		int isDutyCross = 0;
		float dutyWhere = 0;

		float startPhase = 0; // 初始相位 [0,1]

		// PWM 相关参数
		float duty = 0.5f;
		float slope_diff = 8.0f; // 预计算的斜率变化绝对值

		// 辅助内联函数：获取理想无带限的波形值 [-1, 1]
		inline float GetNaiveValue(float p) const
		{
			p -= std::floor(p); // 确保在 [0, 1) 内
			if (p < duty)
				return -1.0f + 2.0f * p / duty; // 上升段
			else
				return 1.0f - 2.0f * (p - duty) / (1.0f - duty); // 下降段
		}

		// 辅助内联函数：获取当前的理论斜率 (每相位单位)
		inline float GetNaiveSlope(float p) const
		{
			p -= std::floor(p);
			if (p < duty) return 2.0f / duty;
			else return -2.0f / (1.0f - duty);
		}

	public:
		TriOscillator()
		{
			SetPWM(0.5f); // 默认完美对称三角波
		}

		// 新增：设置 PWM (占空比)
		inline void SetPWM(float d)
		{
			// 限制在安全范围内，防止斜率为无穷大导致除零崩溃
			if (d < 0.001f) d = 0.001f;
			if (d > 0.999f) d = 0.999f;
			duty = d;

			// 预计算斜率跳变的绝对差值
			// 上升斜率 = 2/duty, 下降斜率 = -2/(1-duty)
			// 差值 = (2/duty) - (-2/(1-duty)) = 2 / (duty * (1 - duty))
			slope_diff = 2.0f / (duty * (1.0f - duty));
		}

		inline void SetStartPhase(float phi)
		{
			startPhase = phi;
		}

		inline void SetPhase(float phi, float where = 0)
		{
			float newt = phi;
			// 强行跳相，需要同时修补幅度 (BLEP) 和 斜率 (BLAMP)
			float val_diff = GetNaiveValue(newt) - GetNaiveValue(t);
			float slope_d = GetNaiveSlope(newt) - GetNaiveSlope(t);

			if (std::abs(val_diff) > 1e-6f) blep.Add(val_diff, where, 1);
			if (std::abs(slope_d) > 1e-6f) blep.Add(slope_d * dt, where, 2);

			t = newt;
		}

		inline bool IsWrapThisSample() const final override { return isWrap; }
		inline float GetWrapWhere() const final override { return where2; }
		inline float GetDT() const final override { return dt; }

		inline void Step(float dt1) final override
		{
			dt = dt1;
			if (dt > 1.0f) dt = 1.0f;
			if (dt < -1.0f) dt = -1.0f;

			float t_old = t;
			t += dt;

			isWrap = 0;
			isDutyCross = 0;

			if (dt > 0.0f)
			{
				if (t >= 1.0f) WrapPhaseDown();
				// 检测是否跨越了顶部峰值 (duty)
				if (t_old < duty && t >= duty)
				{
					isDutyCross = 1;
					dutyWhere = (t - duty) / dt;
				}
			}
			else if (dt < 0.0f)
			{
				if (t < 0.0f) WrapPhaseUp();
				// 反向检测跨越顶部峰值
				if (t_old > duty && t <= duty)
				{
					isDutyCross = 1;
					dutyWhere = (t - duty) / dt;
				}
			}
		}

		inline void WrapPhaseDown()
		{
			t2 = t - 1.0f;
			where2 = t2 / dt;
			if (where2 < 0.0f) where2 = 0.0f;
			if (where2 > 1.0f) where2 = 1.0f;
			isWrap = 1;
		}

		inline void WrapPhaseUp()
		{
			t2 = t + 1.0f;
			where2 = t2 / dt;
			if (where2 < 0.0f) where2 = 0.0f;
			if (where2 > 1.0f) where2 = 1.0f;
			isWrap = 1;
		}

		inline void ApplyWrap()
		{
			t = t2;
			// 底部峰值折角：斜率从下降突变为上升
			// 必须使用 Stage 2 (BLAMP)，且缩放必须乘以 dt！
			blep.Add(slope_diff * dt, where2, 2);
		}

		inline void ApplyDutyCross()
		{
			if (dutyWhere < 0.0f) dutyWhere = 0.0f;
			if (dutyWhere > 1.0f) dutyWhere = 1.0f;
			// 顶部峰值折角：斜率从上升突变为下降 (所以是负号)
			blep.Add(-slope_diff * dt, dutyWhere, 2);
		}

		inline void DoSync(float dstWhere)
		{
			// 1. 如果自己原本的峰值发生在被硬同步之前，先处理自己的
			if (isWrap && dstWhere < where2) {
				ApplyWrap();
				isWrap = 0;
			}
			if (isDutyCross && dstWhere < dutyWhere) {
				ApplyDutyCross();
				isDutyCross = 0;
			}

			// 2. 精确计算同步前后的相位
			float phase_before = (t - dt) + dt * (1.0f - dstWhere);
			float phase_after = startPhase;

			// 3. 极其关键！硬同步会导致波形被“生硬切断”
			// 这意味着既有幅度的突变 (BLEP)，也有斜率的突变 (BLAMP)！
			float val_before = GetNaiveValue(phase_before);
			float val_after = GetNaiveValue(phase_after);
			blep.Add(val_after - val_before, dstWhere, 1); // Stage 1: BLEP

			float slope_before = GetNaiveSlope(phase_before);
			float slope_after = GetNaiveSlope(phase_after);
			blep.Add((slope_after - slope_before) * dt, dstWhere, 2); // Stage 2: BLAMP

			// 4. 计算被同步后，剩余时间跑完的最终相位
			t = phase_after + dstWhere * dt;

			// 5. 检测在同步之后的剩余时间里，是否又跨越了峰值
			if (dt > 0.0f) {
				if (t >= 1.0f) {
					t -= 1.0f;
					blep.Add(slope_diff * dt, (t - 0.0f) / dt, 2);
				}
				else if (phase_after < duty && t >= duty) {
					blep.Add(-slope_diff * dt, (t - duty) / dt, 2);
				}
			}
			else if (dt < 0.0f) {
				if (t < 0.0f) {
					t += 1.0f;
					blep.Add(slope_diff * dt, (t - 1.0f) / dt, 2);
				}
				else if (phase_after > duty && t <= duty) {
					blep.Add(-slope_diff * dt, (t - duty) / dt, 2);
				}
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
			else
			{
				// 如果没有同步打断，常规应用 Wrap 和 顶部峰值的 BLAMP
				if (isWrap) ApplyWrap();
				if (isDutyCross) ApplyDutyCross();
			}

			blep.Step();

			// 因为 GetNaiveValue 直接输出理想的 [-1, 1] 波形
			// 所以直接加上 BLAMP/BLEP 的残差即可
			return GetNaiveValue(t) + blep.Get();
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
			tri = pwm * trifix + tri * 0.995;

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
		TriOscillator osc1;
		TriOscillator osc2;
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
			//osc1.SetWaveform(form);

			osc2.SyncTo(osc1);
			osc2.SetPWM(duty);
			//osc2.SetWaveform(form);

		}
		float ProcessSample()
		{
			osc1.Step(dt1 + fbv * fb);
			osc2.Step(dt2 + fbv * fb);
			float v1 = osc1.Get();
			float v2 = osc2.Get();
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