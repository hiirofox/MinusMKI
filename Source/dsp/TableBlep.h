#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <complex>

namespace TableBlepCoeffs
{
	constexpr static int wsiz = 24;//既决定窗长，又决定阶数
	constexpr static int numTables = 127;
}

class TableBlep
{
private:
	float buf[TableBlepCoeffs::wsiz] = { 0 };
	float v = 0;
	int pos = 0;
	float dc = 0, realdc = 0;
public:
	TableBlep();
	void Add(float amp, float where, int stage = 1);
	void Step();
	float Get();
};