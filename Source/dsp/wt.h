#pragma once

#include <bit>
#include <type_traits> 
#include "IIRBlep.h"
#include "TableBlep.h"

namespace MinusMKI
{
	class TableMutant
	{
	private:
	public:
		virtual void Apply(float* table, int numSamples);
		virtual void SetMutantParams(float param1, float param2, float param3);
	};
	class TableMutantSync :public TableMutant
	{

	};
	class TableMutantSelfm :public TableMutant
	{

	};
	class TableMutantKickizer :public TableMutant
	{

	};
	class TableMutantDisperser :public TableMutant
	{

	};

	class WTOscillator
	{
	public:
		constexpr static int TableWidth = 8192;
	private:
		//IIRBlep2::IIRBlep blit;//只使用其blit功能
		Lagrange4thBlep blit;//其实也不错

		float intMagtable1[TableWidth * 2];
		float intMagtable2[TableWidth * 2];
		float* intMagtable = intMagtable1;
		int swapInterval = 128;
		int sampleCounter = 0;
		int isSwapPrepared = 0;

		float t = 0;
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
		void CalcIntMagtable(float* intMagtable, float* source, int TableWidth)
		{
			int n = ctrz(TableWidth);
			int prevPos = 0;
			int pos = TableWidth;
			int len = TableWidth;
			//const float avge = 1.0 / sqrtf(2.0);
			for (int i = 0; i < TableWidth; ++i)
			{
				intMagtable[i] = source[i];
			}
			for (int m = 0; m < n; ++m)
			{
				int nextLen = len >> 1;
				int curPos = pos;
				//const float avge = 1.0f / cosf(float(M_PI) / float(len));
				float cs = cosf(float(M_PI) / float(len));
				const float avge = 1.0f / (cs * cs);
				for (int i = 0; i < nextLen; ++i)
				{
					/*
					float a = intMagtable[prevPos + i * 2];
					float b = intMagtable[prevPos + i * 2 + 1];
					intMagtable[curPos + i] = (a + b) * avge;
					*/
					int p0 = i * 2;
					int pm = p0 - 1;
					int pp = p0 + 1;
					if (pm < 0)pm += len;
					if (pp >= len)pp -= len;
					float a = intMagtable[prevPos + p0];
					float b = intMagtable[prevPos + pm];
					float c = intMagtable[prevPos + pp];
					intMagtable[curPos + i] = (a + (b + c) * 0.5f) * avge;
				}
				prevPos = curPos;
				pos += nextLen;
				len = nextLen;
			}
			if (pos < TableWidth * 2)
			{
				intMagtable[pos] = intMagtable[prevPos];
			}
		}
		void ProcessSwapTable()
		{
			if (isSwapPrepared == 2)
			{
				intMagtable = (intMagtable == intMagtable1) ? intMagtable2 : intMagtable1;
				isSwapPrepared = 0;
			}
			sampleCounter++;
			if (sampleCounter >= swapInterval)
			{
				sampleCounter = 0;
				isSwapPrepared = 1;
			}
		}
	public:
		WTOscillator()
		{
			float magtable[TableWidth];
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
				magtable[i] = asinf(sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI)) * normv;//tri kick
				//magtable[i] = (sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI) > 0 ? 1.0 : -1.0) * normv;//sqr kick
			}
			CalcIntMagtable(intMagtable, magtable, TableWidth);
		}
		int IsSwapTablePrepared()
		{
			return isSwapPrepared == 1;
		}
		void ApplyTable(float* source, int TableWidth)
		{
			//在单线程调用下一句必死无疑
			//while (!IsSwapTablePrepared())std::this_thread::yield();
			float* nextTable = (intMagtable == intMagtable1) ? intMagtable2 : intMagtable1;
			CalcIntMagtable(nextTable, source, TableWidth);
			isSwapPrepared = 2;
		}
		float ProcessSampleHQ(float dt)//低频模式
		{
			ProcessSwapTable();
			if (dt > 0.499)dt = 0.499;
			if (dt < -0.499)dt = -0.499;

			float absOneDivDt = fabsf(1.0 / dt);

			if (fabsf(dt) < 2.0 / TableWidth)//if(0)
			{
				float posf = t * TableWidth;//TableWidth语境下的t
				int pos1 = posf;
				int pos2 = pos1 + 1;
				float frac = posf - pos1;

				if (pos1 < 0)pos1 += TableWidth;
				if (pos2 < 0)pos2 += TableWidth;

				float mag1 = intMagtable[pos1 % TableWidth];
				float mag2 = intMagtable[pos2 % TableWidth];
				float mag = mag1 + (mag2 - mag1) * frac;

				t += dt;
				if (t >= 1.0)t -= 1.0;
				if (t < 0.0)t += 1.0;
				//blit.Step();
				//return mag * TableWidth + blit.Get();

				blit.Add(mag * TableWidth, 0, 0);
				blit.Step();
				return blit.Get();
			}
			else//高频模式
			{
				int n = ctrz(abs(dt * TableWidth));
				n -= 2; //用前2层表能有效减小跨表相位不连续
				if (n < 0)n = 0;
				float* selectedMagtable = intMagtable + (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
				int selectedTableWidth = TableWidth >> n;

				float ut = t * selectedTableWidth;
				float udt = dt * selectedTableWidth;
				if (udt > 0.0)
				{
					int numpass = (int)(ut + udt) - (int)ut;
					for (int i = 0, pos = ut + 1.0; i < numpass; ++i, ++pos)
					{
						float where = 1.0 - ((float)pos - ut) / udt;
						float mag = selectedMagtable[pos % selectedTableWidth];
						blit.Add(mag * absOneDivDt, where, 0);
					}
				}
				else if (udt < 0.0)
				{
					int numpass = (int)ceilf(ut) - (int)ceilf(ut + udt);
					for (int i = 0, pos = (int)ceilf(ut) - 1; i < numpass; ++i, --pos)
					{
						float where = 1.0 - ((float)pos - ut) / udt;
						int ipos = pos % selectedTableWidth;
						if (ipos < 0)ipos += selectedTableWidth;
						float mag = selectedMagtable[ipos];
						blit.Add(mag * absOneDivDt, where, 0);
					}
				}

				t += dt;
				if (t >= 1.0)t -= 1.0;
				if (t < 0.0)t += 1.0;

				blit.Step();
				return blit.Get();
			}
		}
		float lastt = 0;
		float ProcessSampleLQ(float dt)
		{
			ProcessSwapTable();
			if (dt > 0.499)dt = 0.499;
			if (dt < -0.499)dt = -0.499;

			int n = ctrz(abs(dt * TableWidth));
			float* selectedMagtable = intMagtable + (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
			int selectedTableWidth = TableWidth >> n;
			//float* selectedMagtable = intMagtable; 
			//int selectedTableWidth = TableWidth;

			float posf = t * selectedTableWidth;//TableWidth语境下的t
			int pos1 = posf;
			int pos2 = pos1 + 1;
			float frac = posf - pos1;

			if (pos1 < 0)pos1 += selectedTableWidth;
			if (pos2 < 0)pos2 += selectedTableWidth;

			float mag1 = selectedMagtable[pos1 % selectedTableWidth];
			float mag2 = selectedMagtable[pos2 % selectedTableWidth];
			float mag = mag1 + (mag2 - mag1) * frac;

			lastt = t;
			t += dt;
			if (t >= 1.0)t -= 1.0;
			if (t < 0.0)t += 1.0;
			return mag * selectedTableWidth;
		}
		float ProcessSample(float dt) { return ProcessSampleHQ(dt); }
		void SetStartPhase(float phase) { this->t = phase; }
	};

	class WaveTableOscTest
	{
	private:
		WTOscillator osc1;
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

	class WaveTableOscUnisonTest
	{
	private:
		constexpr static int UnisonNum = 1;
		WTOscillator wav[UnisonNum];
		float ftab[UnisonNum] = { 0 };
		float unitvol = 1.0 / sqrtf(UnisonNum);
	public:
		WaveTableOscUnisonTest()
		{
			for (int i = 0; i < UnisonNum; ++i)
			{
				float randphase = (float)rand() / RAND_MAX;
				wav[i].SetStartPhase(randphase);
			}
		}
		void SetParams(float freq, float detune, float sr = 48000)
		{
			for (int i = 0; i < UnisonNum; ++i)
			{
				float f = freq * (1.0 + ((float)i / UnisonNum - 0.5) * detune * 0.05);
				ftab[i] = f / sr;
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
				sum += wav[i].ProcessSample(ftab[i]);
			}
			sum *= unitvol;
			return sum;
		}
	};
}