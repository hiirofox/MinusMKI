#pragma once

#include <bit>
#include <type_traits> 
#include "IIRBlep.h"
#include "TableBlep.h"

namespace MinusMKI
{
	class WaveTable
	{
	public:
		constexpr static int TableWidth = 8192;
		//constexpr static int TableHeight = 1;
	private:
		//IIRBlep2::IIRBlep blit;//只使用其blit功能
		Lagrange4thBlep blit;//其实也不错
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
			float intn = 0;
			for (int i = 0; i < TableWidth; ++i)
			{
				float x = (float)i / (TableWidth - 1);
				//magtable[i] = (x * 2.0 - 1.0) * normv;//saw
				//magtable[i] = (x < 0.5 ? -1 : 1) * normv;//sqr
				//magtable[i] = asinf(sinf(x * 2.0 * M_PI)) * normv;//tri
				//magtable[i] = (intn += (float)(rand() % 10000) / 10000.0 * (rand() % 2 ? 1 : -1))* 0.1 * normv;
				//magtable[i] = sin(100.0 * powf(x, 0.045) * 2.0 * M_PI) * normv;//sin kick
				//magtable[i] = asinf(sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI)) * normv;//tri kick
				magtable[i] = (sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI) > 0 ? 1.0 : -1.0) * normv;//sqr kick
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
			//dt = -dt;
			if (dt > 0.499)dt = 0.499;
			if (dt < -0.499)dt = -0.499;

			float absOneDivDt = fabsf(1.0 / dt);

			if (fabsf(dt) < 2.0 / TableWidth)
			{
				float posf = t * TableWidth;//TableWidth语境下的t
				int pos1 = posf;
				int pos2 = pos1 + 1;
				float frac = posf - pos1;

				if (pos1 < 0)pos1 += TableWidth;
				if (pos2 < 0)pos2 += TableWidth;

				float mag1 = magtable[pos1 % TableWidth];
				float mag2 = magtable[pos2 % TableWidth];
				float mag = mag1 + (mag2 - mag1) * frac;

				t += dt;
				if (t >= 1.0)t -= 1.0;
				if (t < 0.0)t += 1.0;
				blit.Step();
				return mag * TableWidth + blit.Get();
			}
			else
			{
				int n = ctrz(abs(dt * TableWidth));
				float* selectedMagtable = intMagtable + (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
				int selectedTableWidth = TableWidth >> n;

				float ut = t * selectedTableWidth;//TableWidth语境下的t
				float udt = dt * selectedTableWidth;
				if (udt > 0.0)
				{
					int numpass = (int)(ut + udt) - (int)ut;
					for (int i = 0, pos = ut + 1.0; i < numpass; ++i, ++pos)
					{
						float where = ((float)pos - ut) / udt;
						float mag = selectedMagtable[pos % selectedTableWidth];
						blit.Add(mag * absOneDivDt, 1.0 - where, 0);
					}
				}
				else if (udt < 0.0)
				{
					int numpass = (int)ceilf(ut) - (int)ceilf(ut + udt);
					for (int i = 0, pos = (int)ceilf(ut) - 1; i < numpass; ++i, --pos)
					{
						float where = ((float)pos - ut) / udt;
						int ipos = pos % selectedTableWidth;
						if (ipos < 0)ipos += selectedTableWidth;
						float mag = selectedMagtable[ipos];
						blit.Add(mag * absOneDivDt, 1.0 - where, 0);
					}
				}

				t += dt;
				if (t >= 1.0)t -= 1.0;
				if (t < 0.0)t += 1.0;

				blit.Step();
				return blit.Get();
			}
		}
	};

	class WaveTableOscTest
	{
	private:
		WaveTable osc1;
		float dt = 0, fb = 0;
	public:
		void SetParams(float freq, float fb, float sr = 48000)
		{
			this->fb = fb;
			dt = freq / sr;
		}
		float z = 0.0;
		float ProcessSample()
		{
			return z = osc1.ProcessSample(dt + z * fb);
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