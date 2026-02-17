#pragma once

#include <functional>

/*
template <typename Sample>
class IntegralCalculator
{
private:
	constexpr static Sample accuracy = 1e-4;
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
	constexpr static int tableSize = 4096;
	constexpr static int orders = 2;
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
};*/

#include <vector>
#include <cmath>
#include <functional>
#include <iostream>

class TableADAA
{
private:
	constexpr static int tableSize = 4096; // 增加表大小不会显著增加生成时间
	int orders = 1;                        // 支持动态设置阶数

	// 存储所有阶数的表。
	// tableData[k][i] 表示第 (k+1) 阶原函数在第 i 个点的数值
	std::vector<std::vector<double>> tableData;

	double rangeL = -1.0;
	double rangeR = 1.0;
	double stepSize = 0.0;
	double invStepSize = 0.0;

	// 获取特定阶数、特定位置的插值结果
	double GetIntFunc(double x, int orderIdx)
	{
		// 映射 x 到 index 域
		double idxRaw = (x - rangeL) * invStepSize;

		// 边界处理 (Clamping)
		if (idxRaw < 0.0) idxRaw = 0.0;
		if (idxRaw > tableSize - 1.0000001) idxRaw = tableSize - 1.0000001;

		int pos = static_cast<int>(idxRaw);
		double frac = idxRaw - pos;

		// 获取引用加速访问
		const std::vector<double>& currentTable = tableData[orderIdx];

		// 线性插值 (ADAA本身通常只需要线性插值查找表，因为高阶平滑由差分保证)
		// 如果需要更高精度查找，可在此处改用 Hermite 插值
		return currentTable[pos] + (currentTable[pos + 1] - currentTable[pos]) * frac;
	}

	// 状态历史记录，用于 ADAA 差分计算
	double lx = 0, lfx = 0;
	// 使用 vector 替代固定数组以支持动态阶数，或者使用最大固定大小
	std::vector<double> xtab;
	std::vector<double> ytab;

public:
	TableADAA(int maxOrder = 2) : orders(maxOrder) {
		xtab.resize(orders + 1, 0.0);
		ytab.resize(orders + 1, 0.0);
	}
	TableADAA(std::function<double(double)> nlfunc, double rangeL = -1.0, double rangeR = 1.0, int maxOrder = 2)
	{
		orders = maxOrder;
		xtab.resize(orders + 1, 0.0);
		ytab.resize(orders + 1, 0.0);
		SetNonlinearFunction(nlfunc, rangeL, rangeR);
	}

	// 核心改进：使用 RK4 生成 n 阶原函数表
	void SetNonlinearFunction(std::function<double(double)> func, double rangeL = -1.0, double rangeR = 1.0)
	{
		this->rangeL = rangeL;
		this->rangeR = rangeR;
		this->stepSize = (rangeR - rangeL) / (tableSize - 1);
		this->invStepSize = 1.0 / this->stepSize;

		// 1. 初始化表内存
		tableData.assign(orders, std::vector<double>(tableSize, 0.0));

		// 2. 初始化 RK4 状态向量
		// state[0] 对应 1阶原函数, state[1] 对应 2阶, ...
		// 初始值可以设为0，因为不定积分的常数项在ADAA差分过程中会被抵消
		std::vector<double> state(orders, 0.0);

		// lambda: 定义微分方程组系统 dy/dx = F(x, y)
		// 输入: x, 当前状态向量 y
		// 输出: 导数向量 dydx
		auto SystemDerivative = [&](double x, const std::vector<double>& y) -> std::vector<double> {
			std::vector<double> dydx(orders);
			dydx[0] = func(x); // y1' = f(x)
			for (int i = 1; i < orders; ++i) {
				dydx[i] = y[i - 1]; // y_n' = y_{n-1}
			}
			return dydx;
			};

		// 3. 填充第一个点
		for (int k = 0; k < orders; ++k) tableData[k][0] = state[k];

		// 4. 步进计算 (RK4 Loop)
		double t = rangeL;
		double h = stepSize;

		for (int i = 0; i < tableSize - 1; ++i)
		{
			// --- 标准 RK4 算法 ---
			// k1
			auto k1 = SystemDerivative(t, state);

			// k2
			std::vector<double> state_k2 = state;
			for (int j = 0; j < orders; ++j) state_k2[j] += 0.5 * h * k1[j];
			auto k2 = SystemDerivative(t + 0.5 * h, state_k2);

			// k3
			std::vector<double> state_k3 = state;
			for (int j = 0; j < orders; ++j) state_k3[j] += 0.5 * h * k2[j];
			auto k3 = SystemDerivative(t + 0.5 * h, state_k3);

			// k4
			std::vector<double> state_k4 = state;
			for (int j = 0; j < orders; ++j) state_k4[j] += h * k3[j];
			auto k4 = SystemDerivative(t + h, state_k4);

			// 更新状态
			for (int j = 0; j < orders; ++j) {
				state[j] += (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
				// 存入表
				tableData[j][i + 1] = state[j];
			}

			t += h;
		}
	}

	double Calc(double x)
	{
		// 注意：这是 ADAA 的标准应用逻辑
		// 对于 n 阶 ADAA，我们需要的是第 n 阶原函数的值
		// 然后进行 n 阶有限差分

		// 此处为了兼容你原本的逻辑 (orders <= 1 时直接返回，>1 时做高阶处理)
		// 实际上 ADAA1 需要 1阶表，ADAA2 需要 2阶表
		// 假设 orders 变量代表我们要使用的 ADAA 阶数

		// 获取最高阶的原函数值 (索引是 orders - 1)
		double fx = GetIntFunc(x, orders - 1);

		// 下面这部分是你原来的有限差分逻辑，保持不变
		// 逻辑是：计算斜率 -> 更新历史 -> 如果是高阶则继续计算高阶斜率

		// 注意：ADAA 的标准公式通常是 (AD_n(x) - AD_n(lx)) / (x - lx) 得到 AD_{n-1} 的平均值
		// 你原本的代码逻辑似乎是基于同样的思想。

		double y = (fx - lfx) / (x - lx + 1.0e-18); // 防止除零
		lfx = fx;
		lx = x;

		if (orders <= 1) return y; // ADAA1 完成

		// 移位历史数据
		for (int i = orders; i > 0; --i) {
			xtab[i] = xtab[i - 1];
		}
		xtab[0] = lx;

		// 高阶差分
		// 注意：这里的 y 目前是 (n-1) 阶原函数的近似值
		// 我们需要继续差分直到还原回 0 阶 (原信号)

		// 你的原始逻辑有点独特，标准的 ADAA2 是：
		// out = (AD2(x) - 2*AD2(x-1) + AD2(x-2)) / period^2 ... 类似这种
		// 但针对任意采样率，确实需要除以 delta_x。

		// 让我们沿用你原本的差分逻辑，因为这是正确的“分治差分”（Divided Differences）
		for (int m = 2; m <= orders; ++m)
		{
			double y_prev = ytab[m - 1];
			ytab[m - 1] = y;
			double x_nm = xtab[m];
			double denominator = x - x_nm;
			if (std::abs(denominator) > 1.0e-12)
			{
				// 系数 m 是因为：d(x^n)/dx = n*x^{n-1}
				// 在多次积分后，系数 1/n! 被引入，差分时需要乘回来
				y = (double)m * (y - y_prev) / denominator;
			}
			else {
				// 当 x 变化极小时（例如直流输入），回退到非 ADAA 的计算或者保持上一值
				// 这里为了简单，如果分母太小，往往意味着没有混叠风险，可以直接算原函数
				// 但在 ADAA 类里通常不做这个 fallback，这里仅做保护
				y = y_prev;
			}
		}
		return y;
	}

	// 辅助：重置历史状态（例如在音符触发时）
	void Reset() {
		lfx = 0; lx = 0;
		std::fill(xtab.begin(), xtab.end(), 0.0);
		std::fill(ytab.begin(), ytab.end(), 0.0);
	}
};