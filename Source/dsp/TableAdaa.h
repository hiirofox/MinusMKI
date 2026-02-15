#pragma once

#include <functional>

template <typename Sample>
class IntegralCalculator
{
private:
	constexpr static Sample accuracy = 1e-5;
	std::function<Sample(Sample)> func = [](Sample x) {return x; };
	Sample factorial(int n) {
		Sample res = 1.0;
		for (int i = 2; i <= n; ++i) res *= i;
		return res;
	}
public:
	void SetFunction(std::function<Sample(Sample)> func) { this->func = func; }
	Sample Calc(Sample a, Sample b, int orders = 1,
		int isInCalc = 0, Sample llv = 1e99, Sample av = 0, Sample bv = 0)
	{
		if (isInCalc == 0 && orders > 1)
		{
			auto original_func = func;
			Sample fact = factorial(orders - 1);
			func = [original_func, b, orders, fact](Sample x) {
				return original_func(x) * std::pow(b - x, orders - 1) / fact;
				};
			Sample result = Calc(a, b, 1);
			func = original_func;
			return result;
		}
		Sample mid = (a + b) * 0.5;
		Sample mv = func(mid);
		if (!isInCalc)
		{
			av = func(a);
			bv = func(b);
		}
		Sample intv = (av + 4.0 * mv + bv) * (b - a) / 6.0;
		if (std::abs(intv - llv) <= accuracy) return intv;
		return Calc(a, mid, 1, 1, intv * 0.5, av, mv) +
			Calc(mid, b, 1, 1, intv * 0.5, mv, bv);
	}
};

class TableADAA
{
private:
	constexpr static int tableSize = 32768;
	constexpr static int orders = 3;
	std::function<double(double)> nlfunc = [](double x) {return x; };
	double nltable[tableSize] = { 0 };
	IntegralCalculator<double> calculator;
	double rangeL = -1.0, rangeR = 1.0;
	double GetIntFunc(double x)
	{
		x = x - rangeL;
		x /= (rangeR - rangeL);
		x *= tableSize;
		if (x < 0)x = 0;
		if (x >= tableSize - 1)x = tableSize - 2;
		int pos = x;
		float frac = x - pos;
		return nltable[pos] + (nltable[pos + 1] - nltable[pos]) * frac;
	}
	double xtab[orders + 1] = { 0 };
	double ytab[orders + 1] = { 0 };
public:
	TableADAA() {}
	TableADAA(std::function<double(double)> nlfunc, double rangeL = -1.0, double rangeR = 1.0)
	{
		SetNonlinearFunction(nlfunc, rangeL, rangeR);
	}
	void SetNonlinearFunction(std::function<double(double)> nlfunc, double rangeL = -1.0, double rangeR = 1.0)
	{
		this->nlfunc = nlfunc;
		this->rangeL = rangeL;
		this->rangeR = rangeR;
		calculator.SetFunction(nlfunc);
		double k = rangeR - rangeL;
		double b = rangeL;
		for (int i = 0; i < tableSize; ++i)
		{
			double r = (double)(i + 1) / tableSize;
			r = r * k + b;
			nltable[i] = calculator.Calc(rangeL, r, orders);
		}
	}

	double lx = 0, lfx = 0;
	double Calc(double x)
	{
		double fx = GetIntFunc(x);
		double y = (fx - lfx) / (x - lx);
		lfx = fx;
		lx = x;

		if constexpr (orders <= 1) return y;

		for (int i = orders; i > 0; --i)
		{
			xtab[i] = xtab[i - 1];
		}
		xtab[0] = lx;
		for (int m = 2; m <= orders; ++m)
		{
			double y_prev = ytab[m - 1];
			ytab[m - 1] = y;
			double x_nm = xtab[m];
			double denominator = x - x_nm;
			if (std::abs(denominator) > 1.0e-12)
			{
				y = (double)m * (y - y_prev) / denominator;
			}
		}
		return y;
	}
};