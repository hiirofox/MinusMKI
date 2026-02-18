#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <vector>

#include "TableBlep.h"
#include "TableAdaa.h"


class Blit
{
private:
	float t = 0, dt = 0, sr = 48000;
	int n1 = 1.0, n2 = 1.0;
	const float tinyv = 0.000001;
	float normVolume = 1.0;
	float phi = 0.5;
	inline float dirak(float t, int n1, int n2)
	{
		t *= M_PI;
		float sint = sinf(t);
		if (sint > 0 && sint < tinyv)sint = tinyv;
		if (sint < 0 && sint > -tinyv)sint = -tinyv;
		float v = (sinf(t * n2) - sinf(t * n1)) / sint;
		return v * normVolume;
	}
public:
	void SetParams(float freq, float lowf, float highf, float sr)
	{
		this->dt = freq / sr;
		this->n1 = 1 + (((int)(sr / freq * lowf) >> 1) << 1);
		this->n2 = 1 + (((int)(sr / freq * highf) >> 1) << 1);
		this->sr = sr;
		this->normVolume = freq / sr * 2.0;
	}
	void ProcessBlock(float* outl, float* outr, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float vl = ProcessSample();
			float vr = vl;

			outl[i] = vl;
			outr[i] = vr;
		}
	}
	inline float ProcessSample()
	{
		t += dt;
		t -= (int)t;
		return dirak(t + phi, n1, n2);
	}
	void SetPhi(float phi)
	{
		this->phi = phi;
	}
};

class DisperseOsc
{
private:
	constexpr static int numOscs = 16;
	Blit oscs[numOscs];

public:
	void SetParams(float freq, float curve, float disp, float sr)
	{
		for (int i = 0; i < numOscs; ++i)
		{
			float lf = (float)i / numOscs;
			float hf = (float)(i + 1) / numOscs;
			lf = powf(lf, curve * 5.0 + 1.0);
			hf = powf(hf, curve * 5.0 + 1.0);

			oscs[i].SetParams(freq * (disp * 0.1 * i + 1.0), lf, hf, sr);
		}
	}
	void ProcessBlock(float* outl, float* outr, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float vl = ProcessSample();
			float vr = vl;

			outl[i] = vl;
			outr[i] = vr;
		}
	}

	inline float ProcessSample()
	{
		float v = 0.0;
		for (int i = 0; i < numOscs; ++i)
		{
			v += oscs[i].ProcessSample();
		}
		return v;
	}
	void SetPhi(float phi)
	{
		for (int i = 0; i < numOscs; ++i)
			oscs[i].SetPhi(phi);
	}
};

class WaveformOsc
{
private:
	Blit osc1;
	Blit osc2;
	float intv = 0.0;
	float duty = 0.125;
public:
	void SetParams(float freq, float curve, float disp, float sr)
	{
		osc1.SetParams(freq, curve, disp, sr);
		osc2.SetParams(freq, curve, disp, sr);
		osc1.SetPhi(0);
		osc2.SetPhi(duty);
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
		intv = osc1.ProcessSample() - osc2.ProcessSample() + intv * 0.999;
		return intv;
	}

	void SetDuty(float duty)
	{
		this->duty = duty;
		osc2.SetPhi(duty);
	}
};



class Blep
{
private:
	int wsiz = 2;//单边窗长
	constexpr static int MaxBufLen = 1024;//最多允许64个blep同时运行
	float buf[MaxBufLen] = { 0 };//残差叠加缓冲
	int pos = 0;
	float v = 0;
	float PolyBlep(float t)//经过(0,1),(1,0)的特殊优化过的函数
	{
		if (t < 0 || t > 1) return 0.0;
		float t2 = t * t;
		float t3 = t2 * t;
		float t4 = t2 * t2;
		float t5 = t4 * t;
		float t6 = t3 * t3;
		float t7 = t6 * t;
		float p = -4.36023f * t7
			+ 13.5038f * t6
			- 11.5128f * t5
			- 3.25094f * t4
			+ 7.83529f * t3
			- 0.2393f * t2
			- 2.97566f * t
			+ 0.99986f;
		return p;
	}
	float HermiteBlep(float t)
	{
		if (t < 0.0f || t > 1.0f) return 0.0f;
		t *= 0.5f;
		t += 0.5f;
		float p = t * t * (2.0f * t - 3.0f) + 1.0f;
		return p * 2.0f;
	}
	float LagrangeBlep(float t)
	{
		float x = t;
		if (x >= 2.0f) return 0.0f;
		float x2 = x * x;
		float x3 = x2 * x;
		float x4 = x2 * x2;
		if (x < 1.0f)
		{
			return -0.25f * x4 + 0.6666667f * x3 + 0.5f * x2 - 2.0f * x + 1.0f;
		}
		else
		{
			return 0.0833333f * (x4 - 8.0f * x3 + 22.0f * x2 - 24.0f * x + 8.0f);
		}
	}
	float UsingBlep(float t) { return LagrangeBlep(t); }

	float siBuf[1024];
public:
	void Add(float amp, float where)//amp：发生的阶跃的幅度，where：小数延迟（单位为采样）
	{
		float t = -wsiz + where;//从最左边开始
		for (int i = 0; i < wsiz * 2; ++i)//对整个窗口
		{
			float p;
			if (t < 0) p = UsingBlep(-t) * 0.5 - 1.0;
			else p = -UsingBlep(t) * 0.5;
			buf[(pos + i) % MaxBufLen] += amp * p;//残差
			t += 1.0;
		}
	}
	void Step()
	{
		v = buf[pos];//取值
		buf[pos] = 0;//归零
		pos++;//步进
		if (pos >= MaxBufLen)pos = 0;//维护缓冲区
	}
	float GetBlep()
	{
		return v;//获取blep残差
	}
};

class SiBlep
{
private:
	constexpr static int wsiz = 8;//单边窗长
	float buf[wsiz] = { 0 };//残差叠加缓冲
	float buf2[wsiz] = { 0 };
	int pos = 0;
	float v = 0;
public:
	void Add(float amp, float where)
	{
		float t = -(wsiz >> 1) + where;
		float intv = 0;
		float s = 1.0;
		float const1 = 1.0 / wsiz * 2.0;
		for (int i = 0; i < wsiz; ++i)//对整个窗口
		{
			float w = t * const1;
			float w1 = 1.0 - w * w;
			float wd = w1 * w1;
			s = -s;
			float sc = s / t;
			intv += sc * wd;
			buf2[i] = intv;
			t += 1.0;
		}
		float const2 = 1.0 / intv;
		for (int i = 0, j = pos; i < wsiz; ++i, ++j)//对整个窗口
		{
			float p = buf2[i] * const2 - 1.0;
			if (j >= wsiz)j = 0;
			buf[j] += amp * p;
		}
	}
	void Step()
	{
		v = buf[pos];//取值
		buf[pos] = 0;//归零
		pos++;//步进
		if (pos >= wsiz)pos = 0;//维护缓冲区
	}
	float GetBlep()
	{
		return v;//获取blep残差
	}
};


class BlepTest
{
private:
	TableBlep sb;
	float t = 0;
	float dt = 0;
	float k = 0;
public:
	BlepTest()
	{
		t = (float)(rand() * rand() * rand() * rand() % 10000) / 10000.0;
	}
	void SetParams(float freq, float curve, float disp, float sr)
	{
		dt = freq / sr;
		k = curve * 8 - 4;
		if (dt > 1)dt = 1;
	}
	void ProcessBlock(float* outl, float* outr, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float vl = ProcessSampleSaw();
			float vr = vl;

			outl[i] = vl / 8.0;
			outr[i] = vr / 8.0;
		}
	}

	inline float ProcessSampleSaw()
	{
		t += dt;
		if (t >= 1.0)
		{
			int amp = (int)t;
			float frac = t - amp;
			float where = frac / dt;
			sb.Add(-amp * 2.0, where);
			t = t - amp * 2.0;
		}
		sb.Step();
		float v = sb.Get();
		return  t + v;
	}
	inline float ProcessSampleImp()
	{
		t += dt;
		float imp = 0;
		if (t >= 1.0)
		{
			int amp = (int)t;
			float frac = t - amp;
			float where = frac / dt;
			sb.Add(-amp * 2.0, where, 0);
			t = t - amp * 2.0;
			imp = -amp * 2.0;
		}
		sb.Step();
		float v = sb.Get();
		return  imp + v;
	}
	float tristate = 0;
	float dz = 0;
	inline float ProcessSampleTri()
	{
		if (tristate)
		{
			t += dt;
			if (t >= 1.0)
			{
				float frac = t - 1.0;
				float where = frac / dt;
				sb.Add(-dt * 2.0, where, 2);//1dt -> -1dt
				tristate = !tristate;
				t = 1.0 - frac;
			}
		}
		else
		{
			t -= dt;
			if (t <= -1.0)
			{
				float frac = -t - 1.0;
				float where = frac / dt;
				sb.Add(dt * 2.0, where, 2);
				tristate = !tristate;
				t = -1.0 + frac;
			}
		}
		sb.Step();
		float r = sb.Get();
		float v = t + r;
		float dv = (v - dz) / dt;
		dz = v;
		return v;
	}
};

class AdaaTest
{
private:
	TableADAA adaa{ [](double x) {return x; },-1.0,2.0 };
	double t = 0, dt = 0;
public:
	void SetParams(float freq, float curve, float disp, float sr)
	{
		dt = freq / sr;
	}

	double lastv = 0, lastt = 0;
	inline float ProcessSampleSaw()
	{
		t += dt;

		double t1 = t;
		if (t >= 1.0) t -= (int)t;
		double t2 = t + dt;
		double v1 = 0.5 * t1 * t1;
		double v2 = 0.5 * t2 * t2;
		double dv1 = v2 - v1;
		double dt1 = t2 - t1;
		double y = dv1 / dt1;

		return adaa.Calc(t);
		//return y;
		//return t;
	}
};

class UnisonTest
{
private:
	constexpr static int UnisonNum = 1;
	BlepTest wav[UnisonNum];
	TableADAA adaa{ [](double x) {return atan(x * 10.0); }, -5, 5 ,3 };
	float unitvol = 1.0 / sqrtf(UnisonNum);
public:
	void SetParams(float freq, float curve, float disp, float sr)
	{
		for (int i = 0; i < UnisonNum; ++i)
		{
			float f = freq * (1.0 + ((float)i / UnisonNum - 0.5) * curve * 0.05);
			wav[i].SetParams(f, 0, 0, sr);
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
			sum += wav[i].ProcessSampleSaw();
		}
		sum *= unitvol;
		return sum;
	}
};