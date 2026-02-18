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

class TableBlep
{
private:
	float buf[TableBlepCoeffs::wsiz] = { 0 };
	float v = 0;
	int pos = 0;
public:
	TableBlep();
	void Add(float amp, float where, int stage = 1);
	void Step();
	float Get();
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

		//叠加残差
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
	float Get()
	{
		return v;
	}
};
