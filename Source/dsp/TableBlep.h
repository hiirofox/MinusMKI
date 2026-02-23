#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <complex>

namespace TableBlepCoeffs
{
	constexpr static int wsiz = 24;//既决定窗长，又决定阶数
	constexpr static int numTables = 127;
	constexpr static float bandLimit = 0.95;
}

constexpr static int BLIT_MODE = 0;
constexpr static int BLEP_MODE = 1;
constexpr static int BLAMP_MODE = 2;

class TableBlep
{
private:
	float buf[TableBlepCoeffs::wsiz] = { 0 };
	float v = 0;
	int pos = 0;
public:
	TableBlep();
	void Add(float amp, float where, int stage = BLEP_MODE);
	void Step();
	float Get();
};

class IirDcCompensator
{
private:
	float y1 = 0.0f;
	float y2 = 0.0f;

	float a = 0.999f;//改这个！
	float b = a * 0.999f;
public:
	void Add(float blepDC, float where)
	{
		float s1 = std::powf(a, where);
		float s2 = std::powf(b, where);
		float ds = (s1 / (1.0f - a)) - (s2 / (1.0f - b));
		float dcfix = blepDC / ds;
		y1 += dcfix * s1;
		y2 -= dcfix * s2;
	}
	void Step()
	{
		y1 *= a;
		y2 *= b;
	}
	float Get()
	{
		return y1 + y2;
	}
};

class Lagrange4thBlep
{
private:
	float z1 = 0, z2 = 0, z3 = 0, z4 = 0;
	float v = 0.0;

	inline float Poly4thUnder1(float x)
	{
		x = 1.0f + x * (-2.0f + x * (0.5f + x * (0.6666667f + x * -0.25f)));
		return x * 0.5;
	}

	inline float Poly4thOver1(float x)
	{
		float poly = 8.0f + x * (-24.0f + x * (22.0f + x * (-8.0f + x)));
		return 0.0833333f * poly * 0.5;
	}

	IirDcCompensator dcc;
public:
	void Add(float amp, float where, int stage = 1)
	{
		//生成带限阶跃
		float p1 = -1.0 + Poly4thOver1(2.0f - where);
		float p2 = -1.0 + Poly4thUnder1(1.0f - where);
		float p3 = +1.0 - Poly4thUnder1(where);
		float p4 = +1.0 - Poly4thOver1(where + 1.0f);

		//减去理想阶跃得到残差
		p1 = p1 - 0.0;
		p2 = p2 - 0.0;
		p3 = p3 - 1.0;
		p4 = p4 - 1.0;
		p1 *= amp;
		p2 *= amp;
		p3 *= amp;
		p4 *= amp;
		//叠加残差
		z1 += p1;
		z2 += p2;
		z3 += p3;
		z4 += p4;

		dcc.Add(p1 + p2 + p3 + p4, where);
	}
	void Step()
	{
		v = z1;
		z1 = z2;
		z2 = z3;
		z3 = z4;
		z4 = 0;
		dcc.Step();
	}
	float Get()
	{
		//return v;
		return v - dcc.Get();
	}
};
