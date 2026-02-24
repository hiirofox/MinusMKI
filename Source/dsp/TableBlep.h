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

	float a = 0.995f;//改这个！
	float b = a * 0.999f;
public:
	void Add(float blepDC, float where)
	{
		float s1 = std::powf(a, where);//可以用泰勒展开优化
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

	IirDcCompensator dcc;
public:
	void Add(float amp, float where, int stage = 1)
	{
		float p1 = 0.0, p2 = 0.0, p3 = 0.0, p4 = 0.0;
		if (stage == 1)
		{
			float x = where;
			float x2 = x * x;
			p1 = -1.0f + x2 * (-1.0f / 12.0f + x2 * (1.0f / 24.0f));
			p2 = -25.0f / 24.0f + x2 * (1.0f / 2.0f + x * (1.0f / 6.0f + x * (-1.0f / 8.0f)));
			p3 = -1.0f / 2.0f + x * (1.0f + x * (-1.0f / 4.0f + x * (-1.0f / 3.0f + x * (1.0f / 8.0f))));
			p4 = 1.0f / 24.0f + x2 * (-1.0f / 6.0f + x * (1.0f / 6.0f + x * (-1.0f / 24.0f)));
		}
		else if (stage == 2)
		{
			float x = where;
			p1 = x * (-1.0f + x * x * (41.0f / 144.0f + x * (-15.0f / 128.0f + x * (77.0f / 3840.0f))));
			p2 = -9359.0f / 11520.0f + x * (-395.0f / 768.0f + x * (45.0f / 128.0f + x * (49.0f / 384.0f + x * (-13.0f / 768.0f + x * (-17.0f / 1280.0f)))));
			p3 = -79.0f / 90.0f + x * (7.0f / 16.0f + x * (1.0f / 2.0f + x * (-23.0f / 96.0f + x * (-1.0f / 12.0f + x * (47.0f / 1280.0f)))));
			p4 = -2609.0f / 11520.0f + x * (437.0f / 768.0f + x * (-45.0f / 128.0f + x * (-109.0f / 1152.0f + x * (77.0f / 768.0f + x * (13.0f / 3840.0f)))));
		}

		//叠加残差
		z1 = z1 + p1 * amp;
		z2 = z2 + p2 * amp;
		z3 = z3 + p3 * amp;
		z4 = z4 + p4 * amp;

		//dcc.Add((p1 + p2 + p3 + p4) * amp, where);
	}
	void Step()
	{
		v = z1;
		z1 = z2;
		z2 = z3;
		z3 = z4;
		z4 = 0;
		//dcc.Step();
	}
	float Get()
	{
		return v;
		//return v - dcc.Get();
	}
};
