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

	float a = 0.999f;
	float b = a * 0.999f; 

	float z_a = a - 1.0f;
	float z_a_sq_half = (a - 1.0f) * (a - 1.0f) * 0.5f;
	float z_b = b - 1.0f;
	float z_b_sq_half = (b - 1.0f) * (b - 1.0f) * 0.5f;
	float inv_1_minus_a = 1.0f / (1.0f - a);
	float inv_1_minus_b = 1.0f / (1.0f - b);

public:
	void SetA(float newA)
	{
		a = newA;
		b = a * 0.999f;
		z_a = a - 1.0f;
		z_a_sq_half = z_a * z_a * 0.5f;
		z_b = b - 1.0f;
		z_b_sq_half = z_b * z_b * 0.5f;
		inv_1_minus_a = 1.0f / (1.0f - a);
		inv_1_minus_b = 1.0f / (1.0f - b);
	}

	void Add(float blepDC, float where)
	{
		float x_x_minus_1 = where * (where - 1.0f);
		float s1 = 1.0f + where * z_a + x_x_minus_1 * z_a_sq_half;
		float s2 = 1.0f + where * z_b + x_x_minus_1 * z_b_sq_half;
		float ds = (s1 * inv_1_minus_a) - (s2 * inv_1_minus_b);
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

	inline float SquareWelch(float x)
	{
		x = (x - 2.0f) * 0.5f;
		x = x * x;
		x = 1.0f - x;
		x = x * x;
		return x;
	}
	const float blepIntV = 15.0f / 16.0f;
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
			//p1 += SquareWelch(x + 0.0f) * blepIntV;
			//p2 += SquareWelch(x + 1.0f) * blepIntV;
			//p3 += SquareWelch(x + 2.0f) * blepIntV;
			//p4 += SquareWelch(x + 3.0f) * blepIntV;
			dcc.Add((p1 + p2 + p3 + p4) * amp, where);
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
		return v - dcc.Get();
	}
};
