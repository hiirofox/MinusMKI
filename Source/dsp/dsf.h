#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <vector>

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

class TableBlep
{
private:
	constexpr static int wsiz = 6;//既决定窗长，又决定阶数
	constexpr static int numTables = 12;
	float tableBlit[numTables + 1][wsiz] = { 0 };
	float tableBlep[numTables + 1][wsiz] = { 0 };
	float tableBlamp[numTables + 1][wsiz] = { 0 };
	float buf[wsiz] = { 0 };
	float v = 0;
	int pos = 0;

	void DFT(std::vector<std::complex<double>>& x, int n)
	{
		std::vector<std::complex<double>> y;
		y.resize(n, 0);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				double k = (double)-i * j * 2.0 * M_PI / n;
				std::complex<double> r{ cosf(k),sinf(k) };
				y[i] += x[j] * r;
			}
		}
		for (int i = 0; i < n; ++i)x[i] = y[i];
	}
	void IDFT(std::vector<std::complex<double>>& x, int n)
	{
		std::vector<std::complex<double>> y;
		y.resize(n, 0);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				double k = (double)i * j * 2.0 * M_PI / n;
				std::complex<double> r{ cosf(k),sinf(k) };
				y[i] += x[j] * r;
			}
		}
		for (int i = 0; i < n; ++i)x[i] = y[i] / (double)n;
	}
	float BlackmanHarrisWindow(float x) {
		if (x < 0 || x>1)return 0;
		x = 2.0f * (float)M_PI * x;
		return 0.35875f -
			0.48829f * cosf(x) +
			0.14128f * cosf(2.0f * x) -
			0.01168f * cosf(3.0f * x);
	}

public:
	TableBlep()
	{
		Init(1, 0.75);
	}

	void Init(int usingMinPhase = 1, float wc = 1.0)
	{
		int n = wsiz * (numTables + 1);
		std::vector<std::complex<double>> x;
		x.resize(n, 0);

		for (int i = 0; i < n; ++i)
		{
			float t = (float)i / n * 2.0 - 1.0;
			float wd = BlackmanHarrisWindow((float)i / n);
			float sc;
			float t2 = M_PI * t * wsiz / 2.0 * wc;
			if (fabsf(t2) < 0.000001)sc = 1.0;
			else sc = sin(t2) / (t2);
			x[i] = wd * sc;
		}

		if (usingMinPhase)
		{
			int cepn = n * 4;
			std::vector<std::complex<double>> x2(cepn, 0);
			std::vector<std::complex<double>> cep(cepn, 0);
			for (int i = 0; i < n; ++i) x2[i] = x[i];
			DFT(x2, cepn);
			for (int i = 0; i < cepn; ++i) cep[i] = std::log(std::abs(x2[i]) + 1e-100);
			IDFT(cep, cepn); // 进入倒谱域
			cep[0] = cep[0];//应用因果窗 
			for (int i = 1; i < cepn / 2; ++i) cep[i] *= 2.0;
			cep[cepn / 2] = cep[cepn / 2];
			for (int i = cepn / 2 + 1; i < cepn; ++i) cep[i] = 0.0;
			DFT(cep, cepn); // 变回频域 
			for (int i = 0; i < cepn; ++i) x2[i] = std::exp(cep[i]);
			IDFT(x2, cepn); // 最小相位冲激
			for (int i = 0; i < n; ++i) x[i] = x2[i];
		}

		std::vector<float> mpblit;
		std::vector<float> mpblep;
		std::vector<float> mpblamp;
		mpblit.resize(n, 0);
		mpblep.resize(n, 0);
		mpblamp.resize(n, 0);

		float intv = 0;//积分
		float intv2 = 0;
		for (int i = 0; i < n; ++i)
		{
			intv += x[i].real();
			intv2 += intv;
			mpblit[i] = x[i].real();
			mpblep[i] = intv;
			mpblamp[i] = intv2;
		}
		for (int i = 0; i < n; ++i)
		{
			mpblit[i] = mpblit[i] / intv * (numTables + 1);
			mpblep[i] = mpblep[i] / intv - 1.0;
			mpblamp[i] = (mpblamp[i] / intv2 - (float)i / n) * (numTables + 1);
		}

		for (int i = 0; i < numTables + 1; ++i)
		{
			for (int j = 0; j < wsiz; ++j)
			{
				int k = j * numTables + i;//进行降采样
				tableBlit[i][j] = mpblit[k];
				tableBlep[i][j] = mpblep[k];
				tableBlamp[i][j] = mpblamp[k];
			}
			tableBlit[i][0] -= 1.0;
		}
	}

	void Add(float amp, float where, int stage = 1)//0:blit 1:blep 2:blamp
	{
		float(*table)[wsiz] = tableBlep;
		if (stage == 0)table = tableBlit;
		else if (stage == 2)table = tableBlamp;

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
	float Get()
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
	float k = 0;
public:
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
	inline float ProcessSampleTri()
	{
		if (tristate)
		{
			t += dt;
			if (t > 1.0)
			{
				float frac = t - 1.0;
				float where = frac / dt;
				sb.Add(-dt, where, 2);
				tristate = !tristate;
				t = 1.0 - frac;
			}
		}
		else
		{
			t -= dt;
			if (t < -1.0)
			{
				float frac = -t - 1.0;
				float where = frac / dt;
				sb.Add(dt, where, 2);
				tristate = !tristate;
				t = -1.0 + frac;
			}
		}
		sb.Step();
		float v = sb.Get();
		return  t + v;
	}
};