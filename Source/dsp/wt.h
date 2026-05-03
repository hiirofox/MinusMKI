#pragma once

#include <bit>
#include <type_traits> 
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
		float intMagtable[TableWidth * 2];
		int startPos[TableWidth * 2];
		float t = 0;
	public:
		int ctrz(int x)
		{
			int n = 0;
			while (x > 1)
			{
				x >>= 1;
				++n;
			}
			return n;
		}
		WaveTable()
		{
			const float normv = 1.0 / TableWidth;
			for (int i = 0; i < TableWidth; ++i)
			{
				float x = (float)i / (TableWidth - 1);
				magtable[i] = (x * 2.0 - 1.0) * normv;
			}

			int n = ctrz(TableWidth);
			int prevPos = 0;
			int pos = TableWidth;
			int len = TableWidth;
			//const float avge = 1.0 / sqrtf(2.0);
			for (int i = 0; i < TableWidth; ++i)
			{
				intMagtable[i] = magtable[i];
				startPos[i] = 0;
			}
			for (int m = 0; m < n; ++m)
			{
				int nextLen = len >> 1;
				int curPos = pos;
				const float avge = 1.0f / cosf(float(M_PI) / float(len));
				for (int i = 0; i < nextLen; ++i)
				{
					float a = intMagtable[prevPos + i * 2];
					float b = intMagtable[prevPos + i * 2 + 1];
					intMagtable[curPos + i] = (a + b) * avge;
					startPos[curPos + i] = curPos;
				}
				prevPos = curPos;
				pos += nextLen;
				len = nextLen;
			}
			if (pos < TableWidth * 2)
			{
				intMagtable[pos] = intMagtable[prevPos];
				startPos[pos] = prevPos;
			}
		}
		float ProcessSample(float dt)
		{
			if (dt > 0.5)dt = 0.5;
			int n = ctrz(dt * TableWidth) + 1;
			float* selectedMagtable = intMagtable + (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
			int selectedTableWidth = TableWidth >> n;

			float ut = t * selectedTableWidth;//TableWidth语境下的t
			float udt = dt * selectedTableWidth;
			int numpass = (int)(ut + udt) - (int)ut;
			for (int i = 0, pos = ut + 1.0; i < numpass; ++i, ++pos)
			{
				float where = ((float)pos - ut) / udt;
				float mag = selectedMagtable[pos % selectedTableWidth];
				where = where > 1.0 ? 1.0 : where;
				where = where < 0.0 ? 0.0 : where;
				blit.Add(mag / dt, 1.0 - where, 0);
			}

			t += dt;
			if (t >= 1.0)t -= 1.0;
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