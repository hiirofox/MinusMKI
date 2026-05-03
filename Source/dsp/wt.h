#pragma once

#include "IIRBlep.h"

namespace MinusMKI
{
	class WaveTable
	{
	public:
		constexpr static int TableWidth = 8192;
		//constexpr static int TableHeight = 1;
	private:
		IIRBlep2::IIRBlep blit;//只使用其blit功能
		float magtable[TableWidth];
		float t = 0;
	public:
		WaveTable()
		{
			for (int i = 0; i < TableWidth; ++i)
			{
				float x = (float)i / (TableWidth - 1);
				magtable[i] = x * 2.0 - 1.0;
			}
		}
		float ProcessSample(float dt)
		{
			if (dt > 0.5)dt = 0.5;
			dt *= TableWidth;
			float nextt = t + dt;
			int numpass = (int)nextt - (int)t;
			for (int i = 0, pos = t + 1.0; i < numpass; ++i, ++pos)
			{
				float where = ((float)pos - t) / dt;
				float mag = magtable[pos % TableWidth];
				where = where > 1.0 ? 1.0 : where;
				where = where < 0.0 ? 0.0 : where;
				blit.Add(mag, 1.0 - where, 0);
			}
			t = nextt;
			if (t >= TableWidth)t -= TableWidth;
			blit.Step();
			return blit.Get();
		}
	};
	class WaveTableOscTest
	{
	private:
		WaveTable osc1;
		float dt = 0;
	public:
		void SetParams(float freq, float sr = 48000)
		{
			dt = freq / sr;
		}
		float ProcessSample()
		{
			return osc1.ProcessSample(dt);
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
		}
	};

}