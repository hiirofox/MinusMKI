#pragma once

#include "TableBlep.h"
#include <functional>

namespace MinusMKI
{
	class Oscillator
	{
	private:
	protected:
		Oscillator* syncDst = nullptr;
	public:
		void SyncTo(Oscillator& dst)
		{
			syncDst = &dst;
		}
		void UnregSync()
		{
			syncDst = nullptr;
		}

		virtual bool IsWrapThisSample() = 0;//这个采样是否发生wrap
		virtual float GetWrapWhere() = 0;//获取目前的采样在哪里wrap

		virtual void Step(float dt) = 0;//分开是为了更好处理振荡器内复杂状态
		virtual void UpdateState() = 0;
		virtual float Get() = 0;

	};

	class SawOscillator :public Oscillator
	{
	private:
		TableBlep blep;
		float t = 0;
		float dt = 0;
	public:
		bool IsWrapThisSample() override
		{
			return 0;
		}
		float GetWrapWhere() override
		{
			return -1;
		}
		void Step(float dt1) override
		{
			dt = dt1;
			if (dt > 10.0)dt = 10.0;
			if (dt < -10.0)dt = -10.0;
			t += dt;
		}
		void WrapPhaseDown()
		{
			do {
				t -= 1.0;
				float where = t / dt;
				blep.Add(-1.0, where);
			} while (t >= 1.0);
		}
		void WrapPhaseUp()
		{
			do {
				float where = t / dt;
				blep.Add(1.0, where);
				t += 1.0;
			} while (t < 0.0);
		}
		void UpdateState() override
		{
			if (t >= 1.0)WrapPhaseDown();
			else if (t < 0.0)WrapPhaseUp();
		}
		float Get() override
		{
			blep.Step();
			float v = t + blep.Get();
			return v * 2.0 - 1.0;
		}
	};


	class OscTest
	{
	private:
		SawOscillator osc1;
		SawOscillator osc2;
		float dt1 = 0, dt2 = 0;
		float fb = 0;
		float fbv = 0;
	public:
		void SetParams(float freq, float sync, float fb, float sr)
		{
			dt1 = freq / sr;
			dt2 = freq * sync / sr;
			this->fb = fb;
			osc2.SyncTo(osc1);
		}
		float ProcessSample()
		{
			osc1.Step(dt1 + fbv * fb);
			osc2.Step(dt2 + fbv * fb);
			osc1.UpdateState();
			osc2.UpdateState();
			float v1 = osc1.Get();
			float v2 = osc2.Get();
			fbv = v2;
			return fbv;
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