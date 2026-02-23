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
			Oscillator::UnregSync();
			saw1.UnregSync();
			saw2.UnregSync();
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
	class TriOscillator final : public Oscillator
	{
	private:
		Blep blep;
		float t = 0;
		float dt = 0;
		float t_step_start = 0;

		int isWrap = 0;
		float t2 = 0, where2 = 0;

		int isDutyCross = 0;
		float dutyWhere = 0;

		float startPhase = 0;
		float duty = 0.5f;
		float slope_diff = 8.0f;

		int saw_mode = 0;

		inline float GetNaiveValue(float p) const
		{
			p -= std::floor(p);
			if (saw_mode == -1) return 1.0f - 2.0f * p; // 下降锯齿
			if (saw_mode == 1) return -1.0f + 2.0f * p; // 上升锯齿

			if (p < duty)
				return -1.0f + 2.0f * p / duty;
			else
				return 1.0f - 2.0f * (p - duty) / (1.0f - duty);
		}

		inline float GetNaiveSlope(float p) const
		{
			if (saw_mode == -1) return -2.0f;
			if (saw_mode == 1) return 2.0f;

			p -= std::floor(p);
			if (p < duty) return 2.0f / duty;
			else return -2.0f / (1.0f - duty);
		}

	public:
		TriOscillator()
		{
			SetPWM(0.5f);
		}
		inline void SetPWM(float d)
		{
			if (d <= 0.00001f) {
				duty = 0.0f;
				saw_mode = -1; // 纯下降锯齿
			}
			else if (d >= 0.99999f) {
				duty = 1.0f;
				saw_mode = 1;  // 纯上升锯齿
			}
			else {
				duty = d;
				saw_mode = 0;  // 正常三角波
				slope_diff = 2.0f / (duty * (1.0f - duty));
			}
		}

		inline void SetStartPhase(float phi)
		{
			startPhase = phi;
		}

		inline void SetPhase(float phi, float where = 0)
		{
			float newt = phi;
			float val_diff = GetNaiveValue(newt) - GetNaiveValue(t);
			float slope_d = (GetNaiveSlope(newt) - GetNaiveSlope(t)) * dt;

			if (fabsf(val_diff) > 1e-6f) blep.Add(val_diff, where, 1);
			if (fabsf(slope_d) > 1e-6f) blep.Add(slope_d, where, 2);

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

			t_step_start = t;
			t += dt;

			isWrap = 0;
			isDutyCross = 0;

			if (dt > 0.0f)
			{
				if (t >= 1.0f) {
					isWrap = 1;
					t2 = t - 1.0f;
					where2 = t2 / dt;
				}
				if (saw_mode == 0) {
					if (t_step_start < duty && t >= duty) {
						isDutyCross = 1;
						dutyWhere = (t - duty) / dt;
					}
					else if (t_step_start < 1.0f + duty && t >= 1.0f + duty) {
						isDutyCross = 1;
						dutyWhere = (t - (1.0f + duty)) / dt;
					}
				}
			}
			else if (dt < 0.0f)
			{
				if (t < 0.0f) {
					isWrap = 1;
					t2 = t + 1.0f;
					where2 = t / dt;
				}
				if (saw_mode == 0) {
					if (t_step_start > duty && t <= duty) {
						isDutyCross = 1;
						dutyWhere = (t - duty) / dt;
					}
					else if (t_step_start > duty - 1.0f && t <= duty - 1.0f) {
						isDutyCross = 1;
						dutyWhere = (t - (duty - 1.0f)) / dt;
					}
				}
			}
		}

		inline void ApplyWrap()
		{
			t = t2;
			if (saw_mode == 0) {
				blep.Add(slope_diff * fabsf(dt), where2, 2);
			}
			else if (saw_mode == 1) {
				blep.Add(-2.0f, where2, 1);
			}
			else if (saw_mode == -1) {
				blep.Add(2.0f, where2, 1);
			}
		}

		inline void ApplyDutyCross()
		{
			if (dutyWhere < 0.0f) dutyWhere = 0.0f;
			if (dutyWhere > 1.0f) dutyWhere = 1.0f;
			blep.Add(-slope_diff * fabsf(dt), dutyWhere, 2);
		}

		inline void DoSync(float dstWhere)
		{
			if (isWrap && where2 > dstWhere) {
				ApplyWrap();
				isWrap = 0;
			}
			if (isDutyCross && dutyWhere > dstWhere) {
				ApplyDutyCross();
				isDutyCross = 0;
			}

			float phase_before = t_step_start + dt * (1.0f - dstWhere);
			float phase_after = startPhase;

			float val_before = GetNaiveValue(phase_before);
			float val_after = GetNaiveValue(phase_after);
			if (fabsf(val_after - val_before) > 1e-6f) {
				blep.Add(val_after - val_before, dstWhere, 1);
			}

			float slope_before = GetNaiveSlope(phase_before) * dt;
			float slope_after = GetNaiveSlope(phase_after) * dt;
			if (fabsf(slope_after - slope_before) > 1e-6f) {
				blep.Add(slope_after - slope_before, dstWhere, 2);
			}

			float p_start = phase_after;
			float p_end = phase_after + dstWhere * dt;
			t = p_end;

			if (dt > 0.0f) {
				if (p_end >= 1.0f) {
					t -= 1.0f;
					if (saw_mode == 0) blep.Add(slope_diff * fabsf(dt), (p_end - 1.0f) / dt, 2);
					else if (saw_mode == 1) blep.Add(-2.0f, (p_end - 1.0f) / dt, 1);
					else if (saw_mode == -1) blep.Add(2.0f, (p_end - 1.0f) / dt, 1);
				}
				if (saw_mode == 0) {
					if (p_start < duty && p_end >= duty) {
						blep.Add(-slope_diff * fabsf(dt), (p_end - duty) / dt, 2);
					}
					else if (p_start < 1.0f + duty && p_end >= 1.0f + duty) {
						blep.Add(-slope_diff * fabsf(dt), (p_end - (1.0f + duty)) / dt, 2);
					}
				}
			}
			else if (dt < 0.0f) {
				if (p_end < 0.0f) {
					t += 1.0f;
					if (saw_mode == 0) blep.Add(slope_diff * fabsf(dt), (p_end - 0.0f) / dt, 2);
					else if (saw_mode == 1) blep.Add(-2.0f, (p_end - 0.0f) / dt, 1);
					else if (saw_mode == -1) blep.Add(2.0f, (p_end - 0.0f) / dt, 1);
				}
				if (saw_mode == 0) {
					if (p_start > duty && p_end <= duty) {
						blep.Add(-slope_diff * fabsf(dt), (p_end - duty) / dt, 2);
					}
					else if (p_start > duty - 1.0f && p_end <= duty - 1.0f) {
						blep.Add(-slope_diff * fabsf(dt), (p_end - (duty - 1.0f)) / dt, 2);
					}
				}
			}

			isWrap = 0;
			isDutyCross = 0;
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
				if (isWrap) ApplyWrap();
				if (isDutyCross) ApplyDutyCross();
			}

			blep.Step();
			return GetNaiveValue(t) + blep.Get();
		}
	};

	class WaveformOsc2 final : public Oscillator
	{
	private:
		// 核心引擎：一个完美的带限三角波发生器
		TriOscillator triOsc;

		float duty = 0.25f;
		float dt = 0.0f;
		float form = 0.0f; // imp(0) -> pwm(1) -> tri(2)

		// 求导专用的状态寄存器
		float last_tri = 0.0f;
		float last_pwm = 0.0f;

	public:
		WaveformOsc2()
		{
			UnregSync();
			SetPWM(0.25f);
		}

		inline void SetStartPhase(float sp)
		{
			// 直接透传给底层 TriOscillator
			triOsc.SetStartPhase(sp);
		}

		inline void SetPWM(float d)
		{
			// 保持你原有的映射习惯，并钳位防止除零
			d *= 0.5f;

			this->duty = d;
			triOsc.SetPWM(this->duty);
		}

		inline void SetWaveform(float form)
		{
			this->form = form;
		}

		void SyncTo(Oscillator& dst) final override
		{
			Oscillator::SyncTo(dst);
			triOsc.SyncTo(dst);
		}

		void UnregSync() final override
		{
			Oscillator::UnregSync();
			triOsc.UnregSync();
		}

		inline bool IsWrapThisSample() const final override
		{
			return triOsc.IsWrapThisSample();
		}

		inline float GetWrapWhere() const final override
		{
			return triOsc.GetWrapWhere();
		}

		inline float GetDT() const final override
		{
			return dt;
		}

		inline void Step(float _dt) final override
		{
			dt = _dt;
			if (dt > 1.0f) dt = 1.0f;
			if (dt < -1.0f) dt = -1.0f;
			triOsc.Step(dt);
		}

		float Get() final override
		{
			float tri = triOsc.Get();
			float diff = tri - last_tri;
			last_tri = tri;
			float pwm = 0.0f;
			if (dt != 0.0f) {
				float pwmdutyfix = 1.0;
				if (duty / dt >= 1.0) pwmdutyfix = duty * (1.0f - duty);
				else pwmdutyfix = dt * (1.0 - dt);
				pwm = (diff / dt) * pwmdutyfix;
			}
			float imp = pwm - last_pwm;
			last_pwm = pwm;
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
			return imp * mix1 + pwm * mix2 + tri * mix3;
		}
	};

	class WaveformOsc3 :public Oscillator
	{
	private:
		Blep triblep;
		Blep blep;

		float duty = 0.5;//duty∈[0.0,0.5]
		float dt = 0.001;
		float t = 0;
		int isWrap1 = 0, isWrap2 = 0;
		float where1 = 0, where2 = 0;
		float slope1 = 0, slope2 = 0;

		float GetNaiveTri(float x)//x∈[0,1]
		{
			return x < duty ? x / duty : (1.0 - x) / (1.0 - duty);
		}
		float CrossDetector(float x, float dx, float threshold)
		{
			if (dx > 0.0f)
			{
				float C = std::floor(x - threshold) + threshold;
				if (C > x - dx) return (x - C) / dx;
			}
			else if (dx < 0.0f)
			{
				float C = std::ceil(x - threshold) + threshold;
				if (C < x - dx) return (x - C) / dx;
			}
			return -1.0f;
		}
	public:
		void SetPWM(float duty)
		{
			this->duty = duty * 0.5;
		}
		inline void SetStartPhase(float phi)
		{
		}
		inline void SetPhase(float phi, float where = 0)
		{
		}

		inline bool IsWrapThisSample() const final override
		{
			return isWrap2;
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

			isWrap1 = 0;
			isWrap2 = 0;
			where1 = CrossDetector(t, dt, duty);
			where2 = CrossDetector(t, dt, 1.0);
			if (where1 >= 0.0 && where1 <= 1.0)
			{
				isWrap1 = 1;
				slope1 = -fabsf(dt) / (duty * (1.0 - duty));
			}
			if (where2 >= 0.0 && where2 <= 1.0)
			{
				isWrap2 = 1;
				slope2 = fabsf(dt) / (duty * (1.0 - duty));
			}

			t -= floorf(t);
		}

		inline void ApplyWrap1()
		{
			triblep.Add(slope1, where1, BLAMP_MODE);
		}
		inline void ApplyWrap2()
		{
			triblep.Add(slope2, where2, BLAMP_MODE);
		}
		inline void DoSync(float syncWhere)
		{
			// 1. 判断在 Sync 发生 "之前" 是否已经发生了自然反转
			// where 越大代表在当前 sample 中发生得越早。
			// 这里严格使用 >。如果刚好 ==，说明自然反转和 Sync 同时发生，Sync 在物理上会覆盖自然反转。
			bool wrap1_before = isWrap1 && (where1 > syncWhere);
			bool wrap2_before = isWrap2 && (where2 > syncWhere);

			// isFalling 用于精准记录 Sync 发生那一瞬间，波形是否处于下降沿
			bool isFalling = false;

			// 按照时间先后顺序（where 大的先发生）执行 Sync 前的 Wrap 并更新斜率状态
			if (wrap1_before && wrap2_before)
			{
				if (where1 > where2) {
					ApplyWrap1(); // 先过 Peak
					ApplyWrap2(); // 后过 Valley
					isFalling = false; // 过完 Valley 进入上升沿
				}
				else {
					ApplyWrap2(); // 先过 Valley
					ApplyWrap1(); // 后过 Peak
					isFalling = true; // 过完 Peak 进入下降沿
				}
			}
			else if (wrap1_before)
			{
				ApplyWrap1();
				isFalling = true; // 跨越 duty 后，必然处于下降沿
			}
			else if (wrap2_before)
			{
				ApplyWrap2();
				isFalling = false; // 跨越 1.0 后，必然处于上升沿
			}
			else
			{
				// 如果在 Sync 之前没有任何自然反转，
				// 说明 Sync 瞬间的斜率状态等于这个采样周期一开始的状态。
				float t_start = t + (isWrap2 ? (dt > 0.0f ? 1.0f : -1.0f) : 0.0f) - dt;
				float p_start = t_start - floorf(t_start);
				isFalling = (p_start >= duty);
			}

			// 2. 精确计算 Sync 发生瞬间的相位
			float t_sync = t - syncWhere * dt;
			t_sync -= floorf(t_sync); // 限定在 [0, 1) 区间内

			// 3. 补偿 Sync 瞬间的数值阶跃 (BLEP)
			float diffv = 0.0f - GetNaiveTri(t_sync);
			triblep.Add(diffv, syncWhere, BLEP_MODE);

			// 4. 补偿 Sync 瞬间的斜率突变 (BLAMP)
			// 使用推导出的 isFalling 状态，完美避开 t_sync 处于 0.0 或 1.0 时的边界误判
			if (isFalling)
			{
				triblep.Add(slope2, syncWhere, BLAMP_MODE);
			}

			// 5. 计算 Sync 重置后，在本采样周期内继续累加的相位
			float dt_after = syncWhere * dt;
			float t_new_unwrapped = dt_after;

			// 6. 检查 Sync 之后到采样点结束的这段时间内，是否又触发了 Wrap
			float syncWrapWhere1 = CrossDetector(t_new_unwrapped, dt_after, duty);
			float syncWrapWhere2 = CrossDetector(t_new_unwrapped, dt_after, 1.0f);

			if (syncWrapWhere1 >= 0.0f && syncWrapWhere1 <= 1.0f)
			{
				triblep.Add(slope1, syncWrapWhere1 * syncWhere, BLAMP_MODE);
			}
			if (syncWrapWhere2 >= 0.0f && syncWrapWhere2 <= 1.0f)
			{
				triblep.Add(slope2, syncWrapWhere2 * syncWhere, BLAMP_MODE);
			}

			// 7. 更新 Oscillator 最终在这个周期的相位
			t = t_new_unwrapped - floorf(t_new_unwrapped);
		}
		float Get() final override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncWhere = isSync ? syncDst->GetWrapWhere() : -1.0f;
			if (isSync)
			{
				DoSync(syncWhere);
			}
			else
			{
				if (isWrap1) ApplyWrap1();
				if (isWrap2) ApplyWrap2();
			}

			triblep.Step();
			float tri = GetNaiveTri(t) + triblep.Get();
			return tri * 2.0 - 1.0;
		}
	};

	class OscTest
	{
	private:
		WaveformOsc3 osc1;
		WaveformOsc3 osc2;
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
			osc1.Step(dt1);
			osc2.Step(dt2);
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