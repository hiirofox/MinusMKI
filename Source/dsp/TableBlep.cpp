#include "TableBlep.h"

namespace TableBlepCoeffs
{

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
	void FFTCore(std::vector<std::complex<double>>& x, int n, bool inverse)
	{
		for (int i = 1, j = 0; i < n; i++)
		{
			int bit = n >> 1;
			for (; j & bit; bit >>= 1)
				j ^= bit;
			j ^= bit;
			if (i < j) std::swap(x[i], x[j]);
		}
		for (int len = 2; len <= n; len <<= 1)
		{
			double ang = 2 * M_PI / len * (inverse ? 1 : -1);
			std::complex<double> wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len)
			{
				std::complex<double> w(1);
				for (int j = 0; j < len / 2; j++)
				{
					std::complex<double> u = x[i + j];
					std::complex<double> v = x[i + j + len / 2] * w;
					x[i + j] = u + v;
					x[i + j + len / 2] = u - v;
					w *= wlen;
				}
			}
		}
		if (inverse)
		{
			for (int i = 0; i < n; ++i) x[i] /= n;
		}
	}
	void FFT(std::vector<std::complex<double>>& x, int n)
	{
		FFTCore(x, n, false);
	}
	void IFFT(std::vector<std::complex<double>>& x, int n)
	{
		FFTCore(x, n, true);
	}

	float BlackmanHarrisWindow(float x) {
		if (x < 0 || x>1)return 0;
		x = 2.0f * (float)M_PI * x;
		return 0.35875f -
			0.48829f * cosf(x) +
			0.14128f * cosf(2.0f * x) -
			0.01168f * cosf(3.0f * x);
	}

	void ApplyDCCompensation(std::vector<float>& res, int size)
	{
		if (size <= 1) return;

		// 1. 计算当前残差向量中存在的 DC 总量 (积分)
		double currentDCSum = 0.0;
		for (int i = 0; i < size; ++i) {
			currentDCSum += res[i];
		}

		// 如果 DC 已经极小，直接返回，避免浮点精度问题
		if (std::abs(currentDCSum) < 1e-9) return;

		// 2. 计算 Blackman-Harris 窗在离散点上的总能量 (积分)
		double windowSum = 0.0;
		for (int i = 0; i < size; ++i) {
			// 将索引映射到 0.0 -> 1.0
			float x = (float)i / (float)(size - 1);
			windowSum += BlackmanHarrisWindow(x);
		}

		// 防止除以零 (虽然对于正常窗函数不可能发生)
		if (std::abs(windowSum) < 1e-9) return;

		// 3. 计算缩放系数 k
		// 我们要求: currentDCSum - (k * windowSum) = 0
		// 所以: k = currentDCSum / windowSum
		double scale = currentDCSum / windowSum;

		// 4. 应用补偿
		// res_new[i] = res[i] - (Window[i] * scale)
		for (int i = 0; i < size; ++i) {
			float x = (float)i / (float)(size - 1);
			float compensation = BlackmanHarrisWindow(x) * (float)scale;

			// 核心操作：减去平滑的直流分量
			res[i] -= compensation;
		}
	}

	int isInit = 0;
	float tableBlit[numTables + 1][wsiz] = { 0 };
	float tableBlep[numTables + 1][wsiz] = { 0 };
	float tableBlamp[numTables + 1][wsiz] = { 0 };

	void ApplySincBlep(int usingMinPhase = 1, float wc = 1.0)
	{
		if (isInit)return;
		isInit = 1;

		const int ntable = wsiz * (numTables + 1);
		const int n = 32768;

		std::vector<std::complex<double>> x;
		x.resize(n, 0);

		for (int i = 0; i < n; ++i)
		{
			double t = (double)i / n * 2.0 - 1.0;
			double wd = BlackmanHarrisWindow((double)i / n);
			double sc;
			double t2 = M_PI * t * wsiz / 2.0 * wc;
			if (fabs(t2) < 0.000001)sc = 1.0;
			else sc = sin(t2) / (t2);
			x[i] = wd * sc;
		}

		if (usingMinPhase)
		{
			int cepn = n * 4;
			std::vector<std::complex<double>> x2(cepn, 0);
			std::vector<std::complex<double>> cep(cepn, 0);
			for (int i = 0; i < n; ++i) x2[i] = x[i];
			FFT(x2, cepn);
			for (int i = 0; i < cepn; ++i) cep[i] = std::log(std::abs(x2[i]) + 1e-100);
			IFFT(cep, cepn); // 进入倒谱域
			cep[0] = cep[0];//应用因果窗 
			for (int i = 1; i < cepn / 2; ++i) cep[i] *= 2.0;
			cep[cepn / 2] = cep[cepn / 2];
			for (int i = cepn / 2 + 1; i < cepn; ++i) cep[i] = 0.0;
			FFT(cep, cepn); // 变回频域 
			for (int i = 0; i < cepn; ++i) x2[i] = std::exp(cep[i]);
			IFFT(x2, cepn); // 最小相位冲激
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
			mpblamp[i] = -0.5 * (-mpblamp[i] / intv2 + (float)i / n) * (numTables + 1);
		}
		//todo:为mpblit,mpblep,mpblamp叠加直流补偿窗
		ApplyDCCompensation(mpblit, n);
		ApplyDCCompensation(mpblep, n);
		ApplyDCCompensation(mpblamp, n);

		for (int i = 0; i < numTables + 1; ++i)
		{
			for (int j = 0; j < wsiz; ++j)
			{
				int k = j * numTables + i;//max=ntable
				double k2 = k * (n - 1) / ntable;//下采样
				double frac = k2 - (int)k2;
				int pos = k2;

				tableBlit[i][j] = mpblit[pos] * (1.0 - frac) + mpblit[pos + 1] * frac;
				tableBlep[i][j] = mpblep[pos] * (1.0 - frac) + mpblep[pos + 1] * frac;
				tableBlamp[i][j] = mpblamp[pos] * (1.0 - frac) + mpblamp[pos + 1] * frac;

			}
			tableBlit[i][0] -= 1.0;
		}
	}
}

TableBlep::TableBlep()
{
	TableBlepCoeffs::ApplySincBlep(true, 0.95);
}

void TableBlep::Add(float amp, float where, int stage)//0:blit 1:blep 2:blamp
{
	float(*table)[TableBlepCoeffs::wsiz] = TableBlepCoeffs::tableBlep;
	if (stage == 0)table = TableBlepCoeffs::tableBlit;
	else if (stage == 2)table = TableBlepCoeffs::tableBlamp;

	float inf = where * TableBlepCoeffs::numTables;
	int in1 = (int)inf;
	int in2 = in1 + 1;
	float frac = inf - in1;
	float k1 = (1.0 - frac) * amp;
	float k2 = frac * amp;

	for (int i = 0, j = pos; i < TableBlepCoeffs::wsiz; ++i, ++j)
	{
		float a = table[in1][i];
		float b = table[in2][i];
		float p = a * k1 + b * k2;
		if (j >= TableBlepCoeffs::wsiz) j = 0;
		buf[j] += p;
	}
}
void TableBlep::Step()
{
	v = buf[pos];
	buf[pos] = 0;
	pos++;
	if (pos >= TableBlepCoeffs::wsiz) pos = 0;
}
float TableBlep::Get()
{
	return v;
}