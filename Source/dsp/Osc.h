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
	public:
		bool IsWrapThisSample() override
		{
			return isWrap;
		}
		float GetWrapWhere() override
		{
			return where2;
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
			//if (where2 < 0 || where2 > 1) assert(0);
			if (where2 < 0.0)where2 = 0.0;//不要相信浮点
			if (where2 > 1.0)where2 = 1.0;
			isWrap = 1;
		}
		void WrapPhaseUp()//只预备状态，不更新
		{
			t2 = t;
			where2 = t2 / dt;
			//if (where2 < 0 || where2 > 1) assert(0);
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
			float syncPhase = dstWhere * dt;
			float diff = syncPhase - t;
			blep.Add(diff, dstWhere);
			t = syncPhase;
		}
		float Get() override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncTime = isSync ? syncDst->GetWrapWhere() : -1.0f;
			float selfTime = isWrap ? where2 : -1.0f;

			if (isWrap && !isSync)
			{
				ApplyWrap();
			}
			else if (!isWrap && isSync)
			{
				DoSync(syncTime);
			}
			else if (isWrap && isSync)
			{
				if (syncTime > selfTime)
				{
					DoSync(syncTime);
				}
				else
				{
					ApplyWrap();
					DoSync(syncTime);
				}
			}

			blep.Step();
			float v = t + blep.Get();
			return v * 2.0 - 1.0;
		}
	};

	class ImpulseOsc :public Oscillator
	{
	private:
		TableBlep blep;
		float t = 0;
		float dt = 0;

		int isWrap = 0;
		float t2 = 0, amp2 = 0, where2 = 0;
	public:
		bool IsWrapThisSample() override
		{
			return isWrap;
		}
		float GetWrapWhere() override
		{
			return where2;
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
			//if (where2 < 0 || where2 > 1) assert(0);
			if (where2 < 0.0)where2 = 0.0;//不要相信浮点
			if (where2 > 1.0)where2 = 1.0;
			isWrap = 1;
		}
		void WrapPhaseUp()//只预备状态，不更新
		{
			t2 = t;
			where2 = t2 / dt;
			//if (where2 < 0 || where2 > 1) assert(0);
			if (where2 < 0.0)where2 = 0.0;
			if (where2 > 1.0)where2 = 1.0;
			amp2 = 1.0;
			t2 += 1.0;
			isWrap = 1;
		}
		void ApplyWrap()
		{
			t = t2;
			blep.Add(amp2, where2, 0);
		}
		void DoSync(float dstWhere)
		{
			float syncPhase = dstWhere * dt;
			float diff = syncPhase - t;
			blep.Add(diff, dstWhere, 0);
			t = syncPhase;
		}

		float Get() override
		{
			bool isSync = syncDst && syncDst->IsWrapThisSample();
			float syncTime = isSync ? syncDst->GetWrapWhere() : -1.0f;
			float selfTime = isWrap ? where2 : -1.0f;

			if (isWrap && !isSync)
			{
				ApplyWrap();
			}
			else if (!isWrap && isSync)
			{
				DoSync(syncTime);
			}
			else if (isWrap && isSync)
			{
				if (syncTime > selfTime)
				{
					DoSync(syncTime);
				}
				else
				{
					ApplyWrap();
					DoSync(syncTime);
				}
			}
			blep.Step();
			return blep.Get();
		}
	};

	class WavefromOsc : public Oscillator
	{
	private:
		TableBlep blep1;
		TableBlep blep2;
	public:
		bool IsWrapThisSample() override
		{
		}
		float GetWrapWhere() override
		{
		}

		void Step(float dt1) override
		{
		}

		void WrapPhaseDownT1()//只预备状态，不更新
		{
		}
		void WrapPhaseUpT1()//只预备状态，不更新
		{
		}



		void ApplyWrap()
		{
		}
		void DoSync(float dstWhere)
		{
		}

		float Get() override
		{
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