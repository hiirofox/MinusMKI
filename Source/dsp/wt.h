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
		virtual void Apply(float* table, int numSamples) {};
		virtual void SetMutantParams(float param1, float param2, float param3) {};
	};
	template<int TableWidth>
	class TableMutantSync :public TableMutant
	{
	public:
		constexpr static float MaxFreqX = 16.0 / 8192.0;
	private:
		float descTable[TableWidth] = { 0 };
		float depth = 0, phase = 0, smooth = 0;
		float clampf01(float x) { return x - floorf(x); }
	public:
		void Apply(float* table, int numSamples) override
		{
			float dt = MaxFreqX * depth;
			float t = 0;
			float st = dt * numSamples;
			float smooth1 = expf(-smooth * 12.0);
			for (int i = 0; i < numSamples; ++i)
			{
				t += dt;
				st += smooth1 * (t - st);//cool!
				float idxf = clampf01(st + phase) * numSamples;
				int idx1 = idxf;
				int idx2 = idx1 + 1;
				idx2 = idx2 >= numSamples ? 0 : idx2;
				float frac = idxf - idx1;
				float mag1 = table[idx1];
				float mag2 = table[idx2];
				float mag = mag1 + (mag2 - mag1) * frac;
				descTable[i] = mag;
			}
			for (int i = 0; i < numSamples; ++i)table[i] = descTable[i];
		}
		void SetMutantParams(float depth, float phase, float smooth) override
		{
			this->depth = depth;
			this->phase = phase;
			this->smooth = smooth;
		}
	};
	template<int TableWidth>
	class TableMutantSelfm :public TableMutant
	{

	};
	template<int TableWidth>
	class TableMutantKickizer :public TableMutant
	{
	public:
	private:
		float descTable[TableWidth] = { 0 };
		float depth = 0, tmix = 0, rate = 0;
		float clampf01(float x) { return x - floorf(x); }
	public:
		void Apply(float* table, int numSamples) override
		{
			float depth1 = depth * 200.0;
			float rate1 = expf(-rate * 6.0);
			for (int i = 0; i < numSamples; ++i)
			{
				float x = (float)i / numSamples;
				float t1 = depth1 * powf(x, rate1);
				float t2 = x;
				float t = t2 + (t1 - t2) * tmix;
				float idxf = clampf01(t) * numSamples;
				int idx1 = idxf;
				int idx2 = idx1 + 1;
				idx2 = idx2 >= numSamples ? 0 : idx2;
				float frac = idxf - idx1;
				float mag1 = table[idx1];
				float mag2 = table[idx2];
				float mag = mag1 + (mag2 - mag1) * frac;
				descTable[i] = mag;
			}
			for (int i = 0; i < numSamples; ++i)table[i] = descTable[i];
		}
		void SetMutantParams(float depth, float tmix, float rate) override
		{
			this->depth = depth;
			this->tmix = tmix;
			this->rate = rate;
		}
	};
	template<int TableWidth>
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

		float intMagtable1[TableWidth * 2] = { 0 };
		float intMagtable2[TableWidth * 2] = { 0 };
		float intMagtable3[TableWidth * 2] = { 0 };
		float* intMagtableA = intMagtable1;
		float* intMagtableB = intMagtable2;
		float* nextIntMagtable = intMagtable3;

		float swapInterval = 1.0 / 1024.0;
		float sampleCounter = 0;
		int isSwapPrepared = 0;

		float t = 0;
		static int ctrz(int x)
		{
			int n = 0;
			while (x > 1)
			{
				x >>= 1;
				++n;
			}
			return n;
		}
		void ProcessSwapTable()
		{
			if (isSwapPrepared == 0)
			{
				sampleCounter += swapInterval;
			}
			else if (isSwapPrepared == 2)//新缓冲区已准备好
			{
				if (nextIntMagtable == intMagtable1)//交换缓冲区
				{
					intMagtableA = intMagtable3;
					intMagtableB = intMagtable1;
					nextIntMagtable = intMagtable2;
				}
				else if (nextIntMagtable == intMagtable2)
				{
					intMagtableA = intMagtable1;
					intMagtableB = intMagtable2;
					nextIntMagtable = intMagtable3;
				}
				else if (nextIntMagtable == intMagtable3)
				{
					intMagtableA = intMagtable2;
					intMagtableB = intMagtable3;
					nextIntMagtable = intMagtable1;
				}
				isSwapPrepared = 0;
				sampleCounter = 0.0;
			}
			if (sampleCounter > 1.0)
			{
				sampleCounter = 1.0;
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
			CalcIntMagtable(intMagtable1, magtable, TableWidth);
			CalcIntMagtable(intMagtable2, magtable, TableWidth);
			CalcIntMagtable(intMagtable3, magtable, TableWidth);
		}

		static void CalcIntMagtable(float* intMagtable, float* source, int TableWidth)
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

		int IsSwapTablePrepared()
		{
			return isSwapPrepared == 1;
		}
		void ApplyIntMagtable(float* source, int TableWidth)
		{
			//在单线程调用下一句必死无疑
			//while (!IsSwapTablePrepared())std::this_thread::yield();
			//CalcIntMagtable(nextIntMagtable, source, TableWidth);
			for (int i = 0; i < TableWidth * 2; ++i)
				nextIntMagtable[i] = source[i];
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
				pos1 %= TableWidth;
				pos2 %= TableWidth;

				float mag1A = intMagtableA[pos1];
				float mag2A = intMagtableA[pos2];
				float magA = mag1A + (mag2A - mag1A) * frac;
				float mag1B = intMagtableB[pos1];
				float mag2B = intMagtableB[pos2];
				float magB = mag1B + (mag2B - mag1B) * frac;
				float mag = magA + (magB - magA) * sampleCounter;

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
				int tableStart = (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
				float* selectedTableA = intMagtableA + tableStart;
				float* selectedTableB = intMagtableB + tableStart;
				int selectedTableWidth = TableWidth >> n;

				float ut = t * selectedTableWidth;
				float udt = dt * selectedTableWidth;
				if (udt > 0.0)
				{
					int numpass = (int)(ut + udt) - (int)ut;
					for (int i = 0, pos = ut + 1.0; i < numpass; ++i, ++pos)
					{
						float where = 1.0 - ((float)pos - ut) / udt;
						int idx = pos % selectedTableWidth;
						float magA = selectedTableA[idx];
						float magB = selectedTableB[idx];
						float mag = magA + (magB - magA) * sampleCounter;
						blit.Add(mag * absOneDivDt, where, 0);
					}
				}
				else if (udt < 0.0)
				{
					int numpass = (int)ceilf(ut) - (int)ceilf(ut + udt);
					for (int i = 0, pos = (int)ceilf(ut) - 1; i < numpass; ++i, --pos)
					{
						float where = 1.0 - ((float)pos - ut) / udt;
						int idx = pos % selectedTableWidth;
						if (idx < 0)idx += selectedTableWidth;
						float magA = selectedTableA[idx];
						float magB = selectedTableB[idx];
						float mag = magA + (magB - magA) * sampleCounter;
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
			int tableStart = (n == 0 ? 0 : TableWidth * 2 - (TableWidth >> (n - 1)));
			float* selectedTableA = intMagtableA + tableStart;
			float* selectedTableB = intMagtableB + tableStart;
			int selectedTableWidth = TableWidth >> n;
			//float* selectedMagtable = intMagtable; 
			//int selectedTableWidth = TableWidth;

			float posf = t * selectedTableWidth;//TableWidth语境下的t
			int pos1 = posf;
			int pos2 = pos1 + 1;
			float frac = posf - pos1;

			if (pos1 < 0)pos1 += selectedTableWidth;
			if (pos2 < 0)pos2 += selectedTableWidth;
			pos1 %= selectedTableWidth;
			pos2 %= selectedTableWidth;

			float mag1A = selectedTableA[pos1];
			float mag2A = selectedTableA[pos2];
			float magA = mag1A + (mag2A - mag1A) * frac;
			float mag1B = selectedTableB[pos1];
			float mag2B = selectedTableB[pos2];
			float magB = mag1B + (mag2B - mag1B) * frac;
			float mag = magA + (magB - magA) * sampleCounter;

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
	public:
		constexpr static int TableWidth = WTOscillator::TableWidth;
	private:
		float tableSource[TableWidth * 2] = { 0 };
		float table[TableWidth * 2] = { 0 };
		TableMutantSync<TableWidth> mutantSync;
		TableMutantKickizer<TableWidth> mutantKickizer;
		WTOscillator osc1;
		float dt = 0;
	public:
		WaveTableOscTest()
		{
			const float normv = 1.0 / TableWidth;
			float intn = 0;
			for (int i = 0; i < TableWidth; ++i)
			{
				float x = (float)i / (TableWidth - 1);
				//tableSource[i] = (x * 2.0 - 1.0) * normv;//saw
				//tableSource[i] = (x < 0.5 ? -1 : 1) * normv;//sqr
				tableSource[i] = asinf(sinf(x * 2.0 * M_PI)) * normv;//tri
				//tableSource[i] = (intn += (float)(rand() % 10000) / 10000.0 * (rand() % 2 ? 1 : -1))* 0.1 * normv;
				//tableSource[i] = sin(100.0 * powf(x, 0.045) * 2.0 * M_PI) * normv;//sin kick
				//tableSource[i] = asinf(sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI)) * normv;//tri kick
				//tableSource[i] = (sinf(100.0 * powf(x, 0.045) * 2.0 * M_PI) > 0 ? 1.0 : -1.0) * normv;//sqr kick
			}
		}
		void SetParams(float freq, float p1, float p2, float p3, float p4, float p5, float p6, float sr = 48000)
		{
			dt = freq / sr;
			if (osc1.IsSwapTablePrepared())
			{
				mutantSync.SetMutantParams(p1, p2, p3);
				mutantKickizer.SetMutantParams(p4, p5, p6);
				for (int i = 0; i < TableWidth; ++i)table[i] = tableSource[i];
				mutantSync.Apply(table, TableWidth);
				mutantKickizer.Apply(table, TableWidth);
				WTOscillator::CalcIntMagtable(table, table, TableWidth);
				osc1.ApplyIntMagtable(table, TableWidth);
			}
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