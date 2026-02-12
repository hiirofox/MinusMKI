#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

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
	constexpr static int wsiz = 4;//单边窗长
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



class LagrangeBlep
{
private:
	float z1 = 0, z2 = 0, z3 = 0, z4 = 0;
	float v = 0.0;

	inline float CalcLagrangePoly_Under1(float x)
	{
		return 1.0f + x * (-2.0f + x * (0.5f + x * (0.6666667f + x * -0.25f)));
	}

	inline float CalcLagrangePoly_Over1(float x)
	{
		float poly = 8.0f + x * (-24.0f + x * (22.0f + x * (-8.0f + x)));
		return 0.0833333f * poly;
	}

public:
	void Add(float amp, float where)//amp：发生的阶跃的幅度，where：小数延迟（单位为采样）
	{
		float p1 = CalcLagrangePoly_Over1(2.0f - where);
		float p2 = CalcLagrangePoly_Under1(1.0f - where);
		float p3 = CalcLagrangePoly_Under1(where);
		float p4 = CalcLagrangePoly_Over1(where + 1.0f);
		p1 = p1 * 0.5 - 1.0;
		p2 = p2 * 0.5 - 1.0;
		p3 = -p3 * 0.5;
		p4 = -p4 * 0.5;
		z1 += p1 * amp;
		z2 += p2 * amp;
		z3 += p3 * amp;
		z4 += p4 * amp;
	}
	void Step()
	{
		v = z1;
		z1 = z2;
		z2 = z3;
		z3 = z4;
		z4 = 0;
	}
	float GetBlep()
	{
		return v;
	}
};

class LagrangeBlep2
{
private:
	alignas(16) float z[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	float v = 0.0f;
public:
	void Add(float amp, float where)
	{
		__m128 v_w = _mm_set1_ps(where);
		__m128 v_base = _mm_setr_ps(2.0f, 1.0f, 0.0f, 1.0f);
		__m128 v_sign = _mm_setr_ps(-1.0f, -1.0f, 1.0f, 1.0f);
		__m128 x = _mm_add_ps(v_base, _mm_mul_ps(v_w, v_sign));
		__m128 c4 = _mm_setr_ps(0.08333333f, -0.25f, -0.25f, 0.08333333f);
		__m128 c3 = _mm_setr_ps(-0.6666667f, 0.6666667f, 0.6666667f, -0.6666667f);
		__m128 c2 = _mm_setr_ps(1.8333333f, 0.5f, 0.5f, 1.8333333f);
		__m128 c1 = _mm_set1_ps(-2.0f);
		__m128 c0 = _mm_setr_ps(0.6666667f, 1.0f, 1.0f, 0.6666667f);
		__m128 poly = _mm_add_ps(c3, _mm_mul_ps(x, c4));
		poly = _mm_add_ps(c2, _mm_mul_ps(x, poly));
		poly = _mm_add_ps(c1, _mm_mul_ps(x, poly));
		poly = _mm_add_ps(c0, _mm_mul_ps(x, poly));
		__m128 v_post_scale = _mm_setr_ps(0.5f, 0.5f, -0.5f, -0.5f);
		__m128 v_post_bias = _mm_setr_ps(-1.0f, -1.0f, 0.0f, 0.0f);
		__m128 p_final = _mm_add_ps(_mm_mul_ps(poly, v_post_scale), v_post_bias);
		__m128 v_amp = _mm_set1_ps(amp);
		__m128 v_increment = _mm_mul_ps(p_final, v_amp);
		__m128 v_z = _mm_load_ps(z);
		v_z = _mm_add_ps(v_z, v_increment);
		_mm_store_ps(z, v_z);
	}
	void Step()
	{
		v = z[0];
		__m128 v_z = _mm_load_ps(z);
		__m128 v_shifted = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v_z), 4));
		_mm_store_ps(z, v_shifted);
	}
	float GetBlep()
	{
		return v;
	}
};

#include <complex>
class TableBlep
{
private:
	constexpr static int wsiz = 40;//既决定窗长，又决定阶数
	constexpr static int numTables = 16;
	float table[numTables + 1][wsiz] = { 0 };
	float buf[wsiz] = { 0 };
	float v = 0;
	int pos = 0;

	void DFT(std::vector<std::complex<float>>& x)
	{
		int n = wsiz * numTables;
		std::vector<std::complex<float>> y;
		y.resize(n, 0);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				float k = -i * j * 2.0 * M_PI / n;
				std::complex<float> r{ cosf(k),sinf(k) };
				y[i] += x[j] * r;
			}
		}
		for (int i = 0; i < n; ++i)x[i] = y[i];
	}
	void IDFT(std::vector<std::complex<float>>& x)
	{
		int n = wsiz * numTables;
		std::vector<std::complex<float>> y;
		y.resize(n, 0);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				float k = i * j * 2.0 * M_PI / n;
				std::complex<float> r{ cosf(k),sinf(k) };
				y[i] += x[j] * r;
			}
		}
		for (int i = 0; i < n; ++i)x[i] = y[i] / (float)n;
	}

public:
	TableBlep()
	{
		//UsingSiBlep();
		UsingMinPhaseBlep();
	}

	void UsingMinPhaseBlep()
	{
		int n = wsiz * numTables;
		std::vector<std::complex<float>> x;
		std::vector<std::complex<float>> cep;
		x.resize(n, 0);
		cep.resize(n, 0);
		for (int i = 0; i < n; ++i)
		{
			float t = (float)(i - n / 2) / numTables;
			float w = t / wsiz * 2.0;
			float w1 = 1.0f - w * w;
			float wd = w1 * w1;
			float sc;
			if (std::abs(t) < 1e-5f) {
				sc = 1.0f;
			}
			else {
				sc = sinf(M_PI * t) / (M_PI * t);
			}
			x[i] = wd * sc;
		}
		DFT(x);
		for (int i = 0; i < n; ++i) cep[i] = std::log(std::abs(x[i]) + 1e-10);
		IDFT(cep);
		for (int i = 1; i < n / 2; ++i) {
			cep[i] *= 2.0f;
		}
		for (int i = n / 2; i < n; ++i) {
			cep[i] = 0.0f;
		}
		DFT(cep);
		for (int i = 0; i < n; ++i)x[i] = std::exp(cep[i]);
		IDFT(x);//最小相位冲激响应

		/*
		float intv = 0;//积分
		for (int i = 0; i < n; ++i)
		{
			intv += x[i].real() - x[0].real();
			x[i] = intv;
		}
		for (int i = 0; i < n; ++i)
		{
			x[i] /= intv;//归一化
		}*/

		for (int i = 0; i < numTables; ++i)
		{/*
			for (int j = 0; j < wsiz; ++j)
			{
				int k = j * numTables + i;//进行降采样
				table[i][j] = x[k].real() - 1.0;//冲激积分得阶跃
			}*/
			float intv = 0;
			for (int j = 0; j < wsiz; ++j)
			{
				int k = j * numTables + i;//进行降采样
				intv += x[k].real();
				table[i][j] = intv;//冲激积分得阶跃
			}
			for (int j = 0; j < wsiz; ++j)
			{
				table[i][j] /= intv;
				table[i][j] -= 1.0;
			}
		}

	}
	void UsingSiBlep()
	{
		for (int i = 0; i <= numTables; ++i)
		{
			float where = (float)i / numTables;
			float intv = 0;
			for (int j = 0; j < wsiz; ++j)
			{
				float t = (float)j - wsiz / 2 + where;
				float w = t / wsiz * 2.0;
				float w1 = 1.0f - w * w;
				float wd = w1 * w1;
				float sc;
				if (std::abs(t) < 1e-5f) {
					sc = 1.0f;
				}
				else {
					sc = sinf(M_PI * t) / (M_PI * t);
				}
				intv += wd * sc;
				table[i][j] = intv;
			}
			for (int j = 0; j < wsiz; ++j)
			{
				table[i][j] = (table[i][j] / intv) - 1.0f;
			}
		}
	}
	void Add(float amp, float where)
	{
		float inf = where * numTables;
		int in1 = (int)inf;
		int in2 = in1 + 1;
		float frac = inf - in1;
		float k1 = (1.0 - frac) * amp;
		float k2 = frac * amp;
		for (int i = 0, j = pos; i < wsiz; ++i, ++j)
		{
			float a = table[in1][i];
			float b = table[in2][i];
			float p = a * k1 + b * k2;
			if (j >= wsiz) j = 0;
			buf[j] += p;
		}
	}
	void Step()
	{
		v = buf[pos];
		buf[pos] = 0;
		pos++;
		if (pos >= wsiz) pos = 0;
	}
	float GetBlep()
	{
		return v;
	}
};

class BlepTest
{
private:
	TableBlep sb;
	float t = 0;
	float dt = 0;
public:
	void SetParams(float freq, float curve, float disp, float sr)
	{
		dt = freq / sr;
		if (dt > 1)dt = 1;
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
		float v = sb.GetBlep();
		return  t + v;
	}
};