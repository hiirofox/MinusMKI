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
			else if (t < 0.0)
			{
				float overshoot = t;
				t += 1.0;
				blep.Add(1.0, overshoot / dt);
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
		float t_step_start = 0;

		int isWrap = 0;
		float t2 = 0, where2 = 0;

		int isDutyCross = 0;
		float dutyWhere = 0;

		float startPhase = 0;
		float duty = 0.5f;
		float slope_diff = 8.0f;

		inline float GetNaiveValue(float p) const
		{
			p -= std::floor(p);
			if (p < duty)
				return -1.0f + 2.0f * p / duty;
			else
				return 1.0f - 2.0f * (p - duty) / (1.0f - duty);
		}

		inline float GetNaiveSlope(float p) const
		{
			p -= std::floor(p);
			if (p < duty) return 2.0f / duty;
			else return -2.0f / (1.0f - duty);
		}

	public:
		TriOscillator()
		{
			SetPWM(0.5f);
		}

		float naiveDuty = 0.5, naiveSlopeDiff = 0.0;
		inline void SetPWM(float d)
		{
			naiveDuty = d;
			naiveSlopeDiff = 2.0f / (naiveDuty * (1.0f - naiveDuty));
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

			duty = naiveDuty;
			float dtfix = fabsf(dt);
			if (duty < dtfix) { duty = dtfix; slope_diff = 2.0f / (duty * (1.0f - duty)); }
			else if (1.0 - duty < dtfix) { duty = 1.0 - dtfix; slope_diff = 2.0f / (duty * (1.0f - duty)); }
			else slope_diff = naiveSlopeDiff;

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
				if (t_step_start < duty && t >= duty) {
					isDutyCross = 1;
					dutyWhere = (t - duty) / dt;
				}
				else if (t_step_start < 1.0f + duty && t >= 1.0f + duty) {
					isDutyCross = 1;
					dutyWhere = (t - (1.0f + duty)) / dt;
				}

			}
			else if (dt < 0.0f)
			{
				if (t < 0.0f) {
					isWrap = 1;
					t2 = t + 1.0f;
					where2 = t / dt;
				}
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

		inline void ApplyWrap()
		{
			t = t2;
			blep.Add(slope_diff * fabsf(dt), where2, 2);

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
					blep.Add(slope_diff * fabsf(dt), (p_end - 1.0f) / dt, 2);
				}
				if (p_start < duty && p_end >= duty) {
					blep.Add(-slope_diff * fabsf(dt), (p_end - duty) / dt, 2);
				}
				else if (p_start < 1.0f + duty && p_end >= 1.0f + duty) {
					blep.Add(-slope_diff * fabsf(dt), (p_end - (1.0f + duty)) / dt, 2);
				}

			}
			else if (dt < 0.0f) {
				if (p_end < 0.0f) {
					t += 1.0f;
					blep.Add(slope_diff * fabsf(dt), (p_end - 0.0f) / dt, 2);
				}
				if (p_start > duty && p_end <= duty) {
					blep.Add(-slope_diff * fabsf(dt), (p_end - duty) / dt, 2);
				}
				else if (p_start > duty - 1.0f && p_end <= duty - 1.0f) {
					blep.Add(-slope_diff * fabsf(dt), (p_end - (duty - 1.0f)) / dt, 2);
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

	class WaveformOsc3 final :public Oscillator
	{
	private:
		SawOscillator saw1, saw2;
		TriOscillator tri;

		float duty = 0.25, dt = 0.0;
		float form = 0;//imp(0)->pwm(1)->tri(2)
		float mix1 = 0, mix2 = 0, mix3 = 0;

		float startPhase = 0;
	public:
		WaveformOsc3()
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
			tri.SetStartPhase(p1);
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
			tri.SetPWM(duty);
		}
		inline void SetWaveform(float form)
		{
			this->form = form;
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
		}
		void SyncTo(Oscillator& dst) final override
		{
			Oscillator::SyncTo(dst);
			saw1.SyncTo(dst);
			saw2.SyncTo(dst);
			tri.SyncTo(dst);
		}
		void UnregSync() final override
		{
			Oscillator::UnregSync();
			saw1.UnregSync();
			saw2.UnregSync();
			tri.UnregSync();
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
			if (dt > 0.999)dt = 0.999;
			if (dt < -0.999)dt = -0.999;
			saw1.Step(dt);
			saw2.Step(dt);
			tri.Step(dt);
		}
		float lastpwm = 0, lastv = 0;
		float Get() final override
		{
			float v1 = saw1.Get();
			float v2 = saw2.Get();

			float v3 = -tri.Get();

			float pwm = v1 - v2;

			float sawmix = duty * 50.0;
			if (sawmix < 1.0)
			{
				v3 = v1 * (1.0 - sawmix) + v3 * sawmix;
				float dv = v1 - lastv;
				pwm = dv * (1.0 - sawmix) + pwm * sawmix;
			}
			lastv = v1;

			float imp = pwm - lastpwm;
			lastpwm = pwm;

			float out = imp * mix1 + pwm * mix2 + v3 * mix3;
			return out;
		}
	};

	class WaveformOsc4 final : public Oscillator
	{
	private:
		Blep s1_blep;
		Blep s2_blep;
		Blep t_blep;

		float s1_t = 0;
		int s1_isWrap = 0;
		float s1_t2 = 0, s1_amp2 = 0, s1_where2 = 0;
		float s1_startPhase = 0;

		float s2_t = 0;
		int s2_isWrap = 0;
		float s2_t2 = 0, s2_amp2 = 0, s2_where2 = 0;
		float s2_startPhase = 0;

		float t_t = 0;
		float t_t_step_start = 0;
		int t_isWrap = 0;
		float t_t2 = 0, t_where2 = 0;
		int t_isDutyCross = 0;
		float t_dutyWhere = 0;
		float t_startPhase = 0;
		float t_duty = 0.5f;
		float t_slope_diff = 8.0f;
		float t_naiveDuty = 0.5f, t_naiveSlopeDiff = 0.0f;

		float duty = 0.25f, dt = 0.0f;
		float form = 0;
		float mix1 = 0, mix2 = 0, mix3 = 0;
		float startPhase = 0;
		float lastpwm = 0, lastv = 0;

		inline float TriGetNaiveValue(float p) const
		{
			p -= std::floor(p);
			if (p < t_duty)
				return -1.0f + 2.0f * p / t_duty;
			else
				return 1.0f - 2.0f * (p - t_duty) / (1.0f - t_duty);
		}

		inline float TriGetNaiveSlope(float p) const
		{
			p -= std::floor(p);
			if (p < t_duty) return 2.0f / t_duty;
			else return -2.0f / (1.0f - t_duty);
		}

	public:
		WaveformOsc4()
		{
			SetPWM(0.25f);
		}

		inline void SetStartPhase(float sp)
		{
			startPhase = sp;
			float p1 = startPhase;
			float p2 = startPhase + duty; 
			p1 -= (int)p1;
			p2 -= (int)p2;

			s1_startPhase = p1;
			s2_startPhase = p2;
			t_startPhase = p1;
		}

		inline void SetPWM(float d)
		{
			d *= 0.5f;
			this->duty = d;
			float p1 = startPhase;
			float p2 = startPhase + duty;
			p1 -= (int)p1;
			p2 -= (int)p2;

			s1_startPhase = p1;
			s2_startPhase = p2;

			t_naiveDuty = d;
			t_naiveSlopeDiff = 2.0f / (t_naiveDuty * (1.0f - t_naiveDuty));
		}

		inline void SetWaveform(float form)
		{
			this->form = form;
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
		}

		inline bool IsWrapThisSample() const final override { return s1_isWrap; }
		inline float GetWrapWhere() const final override { return s1_where2; }
		inline float GetDT() const final override { return dt; }

		inline void Step(float _dt) final override
		{
			dt = _dt;
			if (dt > 0.999f) dt = 0.999f;
			if (dt < -0.999f) dt = -0.999f;

			float s1_dt = dt;
			if (s1_dt > 1.0f) s1_dt = 1.0f; 
			if (s1_dt < -1.0f) s1_dt = -1.0f;
			s1_t += s1_dt;
			s1_isWrap = 0;
			if (s1_t >= 1.0f) {
				s1_t2 = s1_t - 1.0f;
				s1_amp2 = -1.0f;
				s1_where2 = s1_t2 / s1_dt;
				if (s1_where2 < 0.0f) s1_where2 = 0.0f;
				if (s1_where2 > 1.0f) s1_where2 = 1.0f;
				s1_isWrap = 1;
			}
			else if (s1_t < 0.0f) {
				s1_t2 = s1_t;
				s1_where2 = s1_t2 / s1_dt;
				if (s1_where2 < 0.0f) s1_where2 = 0.0f;
				if (s1_where2 > 1.0f) s1_where2 = 1.0f;
				s1_amp2 = 1.0f;
				s1_t2 += 1.0f;
				s1_isWrap = 1;
			}

			float s2_dt = dt;
			if (s2_dt > 1.0f) s2_dt = 1.0f;
			if (s2_dt < -1.0f) s2_dt = -1.0f;
			s2_t += s2_dt;
			s2_isWrap = 0;
			if (s2_t >= 1.0f) {
				s2_t2 = s2_t - 1.0f;
				s2_amp2 = -1.0f;
				s2_where2 = s2_t2 / s2_dt;
				if (s2_where2 < 0.0f) s2_where2 = 0.0f;
				if (s2_where2 > 1.0f) s2_where2 = 1.0f;
				s2_isWrap = 1;
			}
			else if (s2_t < 0.0f) {
				s2_t2 = s2_t;
				s2_where2 = s2_t2 / s2_dt;
				if (s2_where2 < 0.0f) s2_where2 = 0.0f;
				if (s2_where2 > 1.0f) s2_where2 = 1.0f;
				s2_amp2 = 1.0f;
				s2_t2 += 1.0f;
				s2_isWrap = 1;
			}

			float t_dt1 = dt;
			if (t_dt1 > 1.0f) t_dt1 = 1.0f;
			if (t_dt1 < -1.0f) t_dt1 = -1.0f;

			t_duty = t_naiveDuty;
			float dtfix = fabsf(t_dt1);
			if (t_duty < dtfix) { t_duty = dtfix; t_slope_diff = 2.0f / (t_duty * (1.0f - t_duty)); }
			else if (1.0 - t_duty < dtfix) { t_duty = 1.0 - dtfix; t_slope_diff = 2.0f / (t_duty * (1.0f - t_duty)); }
			else t_slope_diff = t_naiveSlopeDiff;

			t_t_step_start = t_t;
			t_t += t_dt1;
			t_isWrap = 0;
			t_isDutyCross = 0;

			if (t_dt1 > 0.0f) {
				if (t_t >= 1.0f) {
					t_isWrap = 1;
					t_t2 = t_t - 1.0f;
					t_where2 = t_t2 / t_dt1;
				}
				if (t_t_step_start < t_duty && t_t >= t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (t_t - t_duty) / t_dt1;
				}
				else if (t_t_step_start < 1.0f + t_duty && t_t >= 1.0f + t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (t_t - (1.0f + t_duty)) / t_dt1;
				}
			}
			else if (t_dt1 < 0.0f) {
				if (t_t < 0.0f) {
					t_isWrap = 1;
					t_t2 = t_t + 1.0f;
					t_where2 = t_t / t_dt1;
				}
				if (t_t_step_start > t_duty && t_t <= t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (t_t - t_duty) / t_dt1;
				}
				else if (t_t_step_start > t_duty - 1.0f && t_t <= t_duty - 1.0f) {
					t_isDutyCross = 1;
					t_dutyWhere = (t_t - (t_duty - 1.0f)) / t_dt1;
				}
			}
		}

		float Get() final override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncTime = isSync ? syncDst->GetWrapWhere() : -1.0f;

			if (isSync) {
				if (s1_isWrap && syncTime < s1_where2) {
					s1_t = s1_t2;
					s1_blep.Add(s1_amp2, s1_where2);
					s1_isWrap = 0;
				}
				float syncPhase = syncTime * dt + s1_startPhase;
				float diff = syncPhase - s1_t;
				s1_blep.Add(diff, syncTime);
				s1_t = syncPhase;
				if (s1_t >= 1.0f) {
					int amp = s1_t;
					s1_t -= amp;
					s1_blep.Add(-amp, s1_t / dt);
				}
				else if (s1_t < 0.0f) {
					float overshoot = s1_t;
					s1_t += 1.0f;
					s1_blep.Add(1.0f, overshoot / dt);
				}
			}
			else if (s1_isWrap) {
				s1_t = s1_t2;
				s1_blep.Add(s1_amp2, s1_where2);
			}
			s1_blep.Step();
			float v1 = (s1_t + s1_blep.Get()) * 2.0f - 1.0f;

			if (isSync) {
				if (s2_isWrap && syncTime < s2_where2) {
					s2_t = s2_t2;
					s2_blep.Add(s2_amp2, s2_where2);
					s2_isWrap = 0;
				}
				float syncPhase = syncTime * dt + s2_startPhase;
				float diff = syncPhase - s2_t;
				s2_blep.Add(diff, syncTime);
				s2_t = syncPhase;
				if (s2_t >= 1.0f) {
					int amp = s2_t;
					s2_t -= amp;
					s2_blep.Add(-amp, s2_t / dt);
				}
				else if (s2_t < 0.0f) {
					float overshoot = s2_t;
					s2_t += 1.0f;
					s2_blep.Add(1.0f, overshoot / dt);
				}
			}
			else if (s2_isWrap) {
				s2_t = s2_t2;
				s2_blep.Add(s2_amp2, s2_where2);
			}
			s2_blep.Step();
			float v2 = (s2_t + s2_blep.Get()) * 2.0f - 1.0f;

			if (isSync) {
				if (t_isWrap && t_where2 > syncTime) {
					t_t = t_t2;
					t_blep.Add(t_slope_diff * fabsf(dt), t_where2, 2);
					t_isWrap = 0;
				}
				if (t_isDutyCross && t_dutyWhere > syncTime) {
					float dw = t_dutyWhere;
					if (dw < 0.0f) dw = 0.0f;
					if (dw > 1.0f) dw = 1.0f;
					t_blep.Add(-t_slope_diff * fabsf(dt), dw, 2);
					t_isDutyCross = 0;
				}

				float phase_before = t_t_step_start + dt * (1.0f - syncTime);
				float phase_after = t_startPhase;

				float val_before = TriGetNaiveValue(phase_before);
				float val_after = TriGetNaiveValue(phase_after);
				if (fabsf(val_after - val_before) > 1e-6f) {
					t_blep.Add(val_after - val_before, syncTime, 1);
				}

				float slope_before = TriGetNaiveSlope(phase_before) * dt;
				float slope_after = TriGetNaiveSlope(phase_after) * dt;
				if (fabsf(slope_after - slope_before) > 1e-6f) {
					t_blep.Add(slope_after - slope_before, syncTime, 2);
				}

				float p_start = phase_after;
				float p_end = phase_after + syncTime * dt;
				t_t = p_end;

				if (dt > 0.0f) {
					if (p_end >= 1.0f) {
						t_t -= 1.0f;
						t_blep.Add(t_slope_diff * fabsf(dt), (p_end - 1.0f) / dt, 2);
					}
					if (p_start < t_duty && p_end >= t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - t_duty) / dt, 2);
					}
					else if (p_start < 1.0f + t_duty && p_end >= 1.0f + t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - (1.0f + t_duty)) / dt, 2);
					}
				}
				else if (dt < 0.0f) {
					if (p_end < 0.0f) {
						t_t += 1.0f;
						t_blep.Add(t_slope_diff * fabsf(dt), (p_end - 0.0f) / dt, 2);
					}
					if (p_start > t_duty && p_end <= t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - t_duty) / dt, 2);
					}
					else if (p_start > t_duty - 1.0f && p_end <= t_duty - 1.0f) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - (t_duty - 1.0f)) / dt, 2);
					}
				}

				t_isWrap = 0;
				t_isDutyCross = 0;
			}
			else {
				if (t_isWrap) {
					t_t = t_t2;
					t_blep.Add(t_slope_diff * fabsf(dt), t_where2, 2);
				}
				if (t_isDutyCross) {
					float dw = t_dutyWhere;
					if (dw < 0.0f) dw = 0.0f;
					if (dw > 1.0f) dw = 1.0f;
					t_blep.Add(-t_slope_diff * fabsf(dt), dw, 2);
				}
			}

			t_blep.Step();
			float v3 = -(TriGetNaiveValue(t_t) + t_blep.Get());

			float pwm = v1 - v2;

			float sawmix = duty * 50.0f;
			if (sawmix < 1.0f)
			{
				v3 = v1 * (1.0f - sawmix) + v3 * sawmix;
				float dv = v1 - lastv;
				pwm = dv * (1.0f - sawmix) + pwm * sawmix;
			}
			lastv = v1;

			float imp = pwm - lastpwm;
			lastpwm = pwm;

			float out = imp * mix1 + pwm * mix2 + v3 * mix3;
			return out;
		}
	};

	class WaveformOsc5 final : public Oscillator
	{
	private:
		Blep s1_blep;
		Blep s2_blep;
		Blep t_blep;

		float s1_t = 0;
		float s1_t_step_start = 0; // 替代原版的 t_t_step_start
		int s1_isWrap = 0;
		float s1_t2 = 0, s1_amp2 = 0, s1_where2 = 0;
		float s1_startPhase = 0;

		float s2_t = 0;
		int s2_isWrap = 0;
		float s2_t2 = 0, s2_amp2 = 0, s2_where2 = 0;
		float s2_startPhase = 0;

		// t_t, t_t2, t_where2, t_isWrap, t_startPhase 等完全重复的变量已被移除

		int t_isDutyCross = 0;
		float t_dutyWhere = 0;
		float t_duty = 0.5f;
		float t_slope_diff = 8.0f;
		float t_naiveDuty = 0.5f, t_naiveSlopeDiff = 0.0f;

		float duty = 0.25f, dt = 0.0f;
		float form = 0;
		float mix1 = 0, mix2 = 0, mix3 = 0;
		float startPhase = 0;
		float lastpwm = 0, lastv = 0;

		inline float TriGetNaiveValue(float p) const
		{
			p -= std::floor(p);
			if (p < t_duty)
				return -1.0f + 2.0f * p / t_duty;
			else
				return 1.0f - 2.0f * (p - t_duty) / (1.0f - t_duty);
		}

		inline float TriGetNaiveSlope(float p) const
		{
			p -= std::floor(p);
			if (p < t_duty) return 2.0f / t_duty;
			else return -2.0f / (1.0f - t_duty);
		}

	public:
		WaveformOsc5()
		{
			SetPWM(0.25f);
		}

		inline void SetStartPhase(float sp)
		{
			startPhase = sp;
			float p1 = startPhase;
			float p2 = startPhase + duty;
			p1 -= (int)p1;
			p2 -= (int)p2;

			s1_startPhase = p1;
			s2_startPhase = p2;
		}

		inline void SetPWM(float d)
		{
			d *= 0.5f;
			this->duty = d;
			float p1 = startPhase;
			float p2 = startPhase + duty;
			p1 -= (int)p1;
			p2 -= (int)p2;

			s1_startPhase = p1;
			s2_startPhase = p2;

			t_naiveDuty = d;
			t_naiveSlopeDiff = 2.0f / (t_naiveDuty * (1.0f - t_naiveDuty));
		}

		inline void SetWaveform(float form)
		{
			this->form = form;
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
		}

		inline bool IsWrapThisSample() const final override { return s1_isWrap; }
		inline float GetWrapWhere() const final override { return s1_where2; }
		inline float GetDT() const final override { return dt; }

		inline void Step(float _dt) final override
		{
			// 统一控制 dt
			dt = _dt;
			if (dt > 0.999f) dt = 0.999f;
			if (dt < -0.999f) dt = -0.999f;

			// 统一处理 S1 与 T 共享相位
			s1_t_step_start = s1_t;
			s1_t += dt;
			s1_isWrap = 0;
			if (s1_t >= 1.0f) {
				s1_t2 = s1_t - 1.0f;
				s1_amp2 = -1.0f;
				s1_where2 = s1_t2 / dt;
				if (s1_where2 < 0.0f) s1_where2 = 0.0f;
				if (s1_where2 > 1.0f) s1_where2 = 1.0f;
				s1_isWrap = 1;
			}
			else if (s1_t < 0.0f) {
				s1_t2 = s1_t;
				s1_where2 = s1_t2 / dt;
				if (s1_where2 < 0.0f) s1_where2 = 0.0f;
				if (s1_where2 > 1.0f) s1_where2 = 1.0f;
				s1_amp2 = 1.0f;
				s1_t2 += 1.0f;
				s1_isWrap = 1;
			}

			// S2相位
			s2_t += dt;
			s2_isWrap = 0;
			if (s2_t >= 1.0f) {
				s2_t2 = s2_t - 1.0f;
				s2_amp2 = -1.0f;
				s2_where2 = s2_t2 / dt;
				if (s2_where2 < 0.0f) s2_where2 = 0.0f;
				if (s2_where2 > 1.0f) s2_where2 = 1.0f;
				s2_isWrap = 1;
			}
			else if (s2_t < 0.0f) {
				s2_t2 = s2_t;
				s2_where2 = s2_t2 / dt;
				if (s2_where2 < 0.0f) s2_where2 = 0.0f;
				if (s2_where2 > 1.0f) s2_where2 = 1.0f;
				s2_amp2 = 1.0f;
				s2_t2 += 1.0f;
				s2_isWrap = 1;
			}

			// T的Duty交叉计算直接使用 s1_t
			t_duty = t_naiveDuty;
			float dtfix = fabsf(dt);
			if (t_duty < dtfix) { t_duty = dtfix; t_slope_diff = 2.0f / (t_duty * (1.0f - t_duty)); }
			else if (1.0f - t_duty < dtfix) { t_duty = 1.0f - dtfix; t_slope_diff = 2.0f / (t_duty * (1.0f - t_duty)); }
			else t_slope_diff = t_naiveSlopeDiff;

			t_isDutyCross = 0;
			if (dt > 0.0f) {
				if (s1_t_step_start < t_duty && s1_t >= t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (s1_t - t_duty) / dt;
				}
				else if (s1_t_step_start < 1.0f + t_duty && s1_t >= 1.0f + t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (s1_t - (1.0f + t_duty)) / dt;
				}
			}
			else if (dt < 0.0f) {
				if (s1_t_step_start > t_duty && s1_t <= t_duty) {
					t_isDutyCross = 1;
					t_dutyWhere = (s1_t - t_duty) / dt;
				}
				else if (s1_t_step_start > t_duty - 1.0f && s1_t <= t_duty - 1.0f) {
					t_isDutyCross = 1;
					t_dutyWhere = (s1_t - (t_duty - 1.0f)) / dt;
				}
			}
		}

		float Get() final override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncTime = isSync ? syncDst->GetWrapWhere() : -1.0f;

			// 【关键】提前截取给T用的状态，避免下面 s1块修改状态导致判断失败
			int shared_isWrap = s1_isWrap;
			float shared_where2 = s1_where2;

			if (isSync) {
				if (s1_isWrap && syncTime < s1_where2) {
					s1_t = s1_t2;
					s1_blep.Add(s1_amp2, s1_where2);
					s1_isWrap = 0;
				}
				float syncPhase = syncTime * dt + s1_startPhase;
				float diff = syncPhase - s1_t;
				s1_blep.Add(diff, syncTime);
				s1_t = syncPhase;
				if (s1_t >= 1.0f) {
					int amp = s1_t;
					s1_t -= amp;
					s1_blep.Add(-amp, s1_t / dt);
				}
				else if (s1_t < 0.0f) {
					float overshoot = s1_t;
					s1_t += 1.0f;
					s1_blep.Add(1.0f, overshoot / dt);
				}
			}
			else if (s1_isWrap) {
				s1_t = s1_t2;
				s1_blep.Add(s1_amp2, s1_where2);
			}
			s1_blep.Step();
			float v1 = (s1_t + s1_blep.Get()) * 2.0f - 1.0f;

			if (isSync) {
				if (s2_isWrap && syncTime < s2_where2) {
					s2_t = s2_t2;
					s2_blep.Add(s2_amp2, s2_where2);
					s2_isWrap = 0;
				}
				float syncPhase = syncTime * dt + s2_startPhase;
				float diff = syncPhase - s2_t;
				s2_blep.Add(diff, syncTime);
				s2_t = syncPhase;
				if (s2_t >= 1.0f) {
					int amp = s2_t;
					s2_t -= amp;
					s2_blep.Add(-amp, s2_t / dt);
				}
				else if (s2_t < 0.0f) {
					float overshoot = s2_t;
					s2_t += 1.0f;
					s2_blep.Add(1.0f, overshoot / dt);
				}
			}
			else if (s2_isWrap) {
				s2_t = s2_t2;
				s2_blep.Add(s2_amp2, s2_where2);
			}
			s2_blep.Step();
			float v2 = (s2_t + s2_blep.Get()) * 2.0f - 1.0f;

			if (isSync) {
				// 使用上面缓存的 share 状态
				if (shared_isWrap && shared_where2 > syncTime) {
					t_blep.Add(t_slope_diff * fabsf(dt), shared_where2, 2);
				}
				if (t_isDutyCross && t_dutyWhere > syncTime) {
					float dw = t_dutyWhere;
					if (dw < 0.0f) dw = 0.0f;
					if (dw > 1.0f) dw = 1.0f;
					t_blep.Add(-t_slope_diff * fabsf(dt), dw, 2);
					t_isDutyCross = 0;
				}

				float phase_before = s1_t_step_start + dt * (1.0f - syncTime);
				float phase_after = s1_startPhase;

				float val_before = TriGetNaiveValue(phase_before);
				float val_after = TriGetNaiveValue(phase_after);
				if (fabsf(val_after - val_before) > 1e-6f) {
					t_blep.Add(val_after - val_before, syncTime, 1);
				}

				float slope_before = TriGetNaiveSlope(phase_before) * dt;
				float slope_after = TriGetNaiveSlope(phase_after) * dt;
				if (fabsf(slope_after - slope_before) > 1e-6f) {
					t_blep.Add(slope_after - slope_before, syncTime, 2);
				}

				// s1块已经修改过了s1_t本身，这里T只需借用p_end产生对应的Blep计算，不再需要手动加减s1_t/t_t
				float p_start = phase_after;
				float p_end = phase_after + syncTime * dt;

				if (dt > 0.0f) {
					if (p_end >= 1.0f) {
						t_blep.Add(t_slope_diff * fabsf(dt), (p_end - 1.0f) / dt, 2);
					}
					if (p_start < t_duty && p_end >= t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - t_duty) / dt, 2);
					}
					else if (p_start < 1.0f + t_duty && p_end >= 1.0f + t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - (1.0f + t_duty)) / dt, 2);
					}
				}
				else if (dt < 0.0f) {
					if (p_end < 0.0f) {
						t_blep.Add(t_slope_diff * fabsf(dt), (p_end - 0.0f) / dt, 2);
					}
					if (p_start > t_duty && p_end <= t_duty) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - t_duty) / dt, 2);
					}
					else if (p_start > t_duty - 1.0f && p_end <= t_duty - 1.0f) {
						t_blep.Add(-t_slope_diff * fabsf(dt), (p_end - (t_duty - 1.0f)) / dt, 2);
					}
				}

				t_isDutyCross = 0;
			}
			else {
				if (shared_isWrap) {
					t_blep.Add(t_slope_diff * fabsf(dt), shared_where2, 2);
				}
				if (t_isDutyCross) {
					float dw = t_dutyWhere;
					if (dw < 0.0f) dw = 0.0f;
					if (dw > 1.0f) dw = 1.0f;
					t_blep.Add(-t_slope_diff * fabsf(dt), dw, 2);
				}
			}

			t_blep.Step();
			// 最终直接共用计算好的 s1_t 即可
			float v3 = -(TriGetNaiveValue(s1_t) + t_blep.Get());

			float pwm = v1 - v2;

			float sawmix = duty * 50.0f;
			if (sawmix < 1.0f)
			{
				v3 = v1 * (1.0f - sawmix) + v3 * sawmix;
				float dv = v1 - lastv;
				pwm = dv * (1.0f - sawmix) + pwm * sawmix;
			}
			lastv = v1;

			float imp = pwm - lastpwm;
			lastpwm = pwm;

			float out = imp * mix1 + pwm * mix2 + v3 * mix3;
			return out;
		}
	};
	class OscTest
	{
	private:
		WaveformOsc5 osc1;
		WaveformOsc5 osc2;
		float dt1 = 0, dt2 = 0;
		float duty = 0.5;
		float fb = 0, fbv = 0;
	public:
		void SetParams(float freq, float sync, float pwm, float form, float fb, float sr)
		{
			dt1 = freq / sr;
			dt2 = freq * sync / sr;
			this->fb = fb * 0.01;
			this->duty = pwm;
			osc1.SetPWM(duty);
			osc1.SetWaveform(form);

			osc2.SyncTo(osc1);
			osc2.SetPWM(duty);
			osc2.SetWaveform(form);
		}
		float ProcessSample()
		{
			osc1.Step(-dt1 + fb * fbv);
			osc2.Step(-dt2 + fb * fbv);
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
		constexpr static int UnisonNum = 50;
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