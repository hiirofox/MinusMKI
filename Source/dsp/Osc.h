#pragma once

#include "TableBlep.h"
#include <functional>

namespace MinusMKI
{
	class Oscillator
	{
	private:
		std::vector<Oscillator*> syncfunc;
	public:
		void RegisterSyncOscillator(Oscillator& osc)
		{
			if (&osc != this) syncfunc.push_back(&osc);
		}
		void UnregisterSyncOscillator(Oscillator& osc)
		{
			syncfunc.erase(std::remove_if(syncfunc.begin(), syncfunc.end(),
				[&osc](const Oscillator* ptr) { return ptr == &osc; }),
				syncfunc.end());
		}
		void Sync(float where)
		{
			for (auto& func : syncfunc)
			{
				func->ResetPhase(where);
			}
		}
		virtual void PrepareSample(float dt) = 0;//分开是为了更好处理振荡器内复杂状态
		virtual float ProcessSample() = 0;

		virtual void ResetPhase(float where) = 0;
	};

	class SawOscillator :public Oscillator
	{
	private:
		TableBlep blep;
		float t = 0, dt = 0;
		float where = 0, syncwhere = 0;
		int isReset = 0, isSyncReset = 0;
	public:
		void PrepareSample(float dt1) override
		{
			dt = dt1;
			if (dt > 1)dt = 1;

			t += dt;
			if (t >= 1.0)
			{
				where = (t - 1) / dt;
				Sync(where);
				isReset = 1;
			}
		}
		void ResetPhase(float where) override
		{
			syncwhere = where;
			isSyncReset = 1;
		}

		void Wrap1()
		{
			blep.Add(-1.0f, where);
			t -= 1.0f;
		}
		void Reset(float rswhere)
		{
			float t0 = rswhere * dt;
			blep.Add(t0 - t, rswhere);
			t = t0;
		}
		float ProcessSample() override
		{
			if (isSyncReset)
			{
				if (isReset && (where > syncwhere))
				{
					Wrap1();
					Reset(syncwhere);
				}
				else
				{
					Reset(syncwhere);
				}
			}
			else if (isReset)
			{
				Wrap1();
			}
			
			isSyncReset = 0;
			isReset = 0;

			blep.Step();
			float v = t + blep.Get();
			return v * 2.0f - 1.0f;
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
			osc1.RegisterSyncOscillator(osc2);
		}
		float ProcessSample()
		{
			osc1.PrepareSample(dt1 + fbv * fb);
			osc2.PrepareSample(dt2 + fbv * fb);
			float v1 = osc1.ProcessSample();
			float v2 = osc2.ProcessSample();
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