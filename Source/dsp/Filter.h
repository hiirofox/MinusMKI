#pragma once
#include <math.h>

namespace MinusMKI
{
	class Filter
	{
	protected:
		inline float cheapCosPi(float x)
		{
			int n = floorf(x);
			float f = x - n;
			float y = 1.0f + f * f * (4.0f * f - 6.0f);
			return (n & 1) ? -y : y;
		}

		inline float cheapSinPi(float x)
		{
			int n = floorf(x);
			float f = x - n;
			float t = f * (1.0f - f);
			float y = t * (3.14159265f + 3.40185714f * t);
			return (n & 1) ? -y : y;
		}
	public:
		virtual void SetSampleRate(float sampleRate) {};
		virtual void SetFilterParams(float cutoff, float reso, float morph) {};
		virtual void Reset() {};
		virtual float ProcessSample(float x) { return 0; };
	};
	class SVFilter12dB :public Filter
	{
	private:
		float sampleRate = 48000.0;
		float z1 = 0, z2 = 0;
		float d0 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;
	public:
		inline float ProcessSample(float in) override
		{
			float x = in - z1 - z2;
			float out = d0 * x + d1 * z1 + d2 * z2;
			z2 += c2 * z1;
			z1 += c1 * x;
			return out;
		}
		void Reset() override { z1 = z2 = 0; }
		void DesignBasicFilter(float cutoff, float reso, float lp2bp2hp)
		{
			constexpr float pi = 3.14159265358979323846f;

			if (cutoff < 40.0)reso = 0.707;
			cutoff /= sampleRate;
			if (cutoff < 0.00005)cutoff = 0.00005;
			if (cutoff > 0.45)cutoff = 0.45;

			//float w0 = 2.0f * pi * cutoff;
			//float cs = cosf(w0);
			//float sn = sinf(w0);
			float cs = cheapCosPi(2.0 * cutoff);
			float sn = cheapSinPi(2.0 * cutoff);

			float alpha = sn / (2.0f * reso); // reso = Q
			float a0 = 1.0f + alpha;
			float b0 = ((1.0f - cs) * 0.5f) / a0;
			float a1 = (-2.0f * cs) / a0;
			float a2 = (1.0f - alpha) / a0;
			c1 = a1 + 2.0f;
			c2 = (1.0f + a1 + a2) / c1;
			float bp_b0 = (sn * 0.5f) / a0;
			float bp_d0 = bp_b0;
			float bp_d1 = (2.0f * bp_b0) / c1;
			float bp_d2 = 0.0f;
			float t = lp2bp2hp * 2.0f;
			if (t < 1.0f)// LP -> BP
			{
				float lp_d0 = ((1.0f - cs) * 0.5f) / a0;
				float lp_d1 = c2;
				float lp_d2 = 1.0f;
				d0 = lp_d0 + t * (bp_d0 - lp_d0);
				d1 = lp_d1 + t * (bp_d1 - lp_d1);
				d2 = lp_d2 + t * (bp_d2 - lp_d2);
			}
			else// BP -> HP
			{
				t -= 1.0f;
				float hp_d0 = ((1.0f + cs) * 0.5f) / a0;
				float hp_d1 = 0.0f;
				float hp_d2 = 0.0f;
				d0 = bp_d0 + t * (hp_d0 - bp_d0);
				d1 = bp_d1 + t * (hp_d1 - bp_d1);
				d2 = bp_d2 + t * (hp_d2 - bp_d2);
			}
		}
		void SetFilterParams(float cutoff, float reso, float morph) override
		{
			DesignBasicFilter(cutoff, reso, morph);
		}
		void SetSampleRate(float sr) override
		{
			sampleRate = sr;
		}
	};
	class SVFilter24dB :public Filter
	{
	private:
		float sampleRate = 48000.0;
		float z1 = 0, z2 = 0, z3 = 0, z4 = 0;
		float d0 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;
	public:
		inline float ProcessSample(float in) override
		{
			float x1 = in - z1 - z2;
			float out1 = d0 * x1 + d1 * z1 + d2 * z2;
			float x2 = out1 - z3 - z4;
			float out2 = d0 * x2 + d1 * z3 + d2 * z4;
			z2 += c2 * z1;
			z1 += c1 * x1;
			z4 += c2 * z3;
			z3 += c1 * x2;
			return out2;
		}
		void Reset() override { z1 = z2 = z3 = z4 = 0; }
		void DesignBasicFilter(float cutoff, float reso, float lp2bp2hp)
		{
			constexpr float pi = 3.14159265358979323846f;

			if (cutoff < 40.0)reso = 0.707;
			cutoff /= sampleRate;
			if (cutoff < 0.00005)cutoff = 0.00005;
			if (cutoff > 0.45)cutoff = 0.45;

			reso = sqrtf(reso);//!

			//float w0 = 2.0f * pi * cutoff;
			//float cs = cosf(w0);
			//float sn = sinf(w0);
			float cs = cheapCosPi(2.0 * cutoff);
			float sn = cheapSinPi(2.0 * cutoff);

			float alpha = sn / (2.0f * reso); // reso = Q
			float a0 = 1.0f + alpha;
			float b0 = ((1.0f - cs) * 0.5f) / a0;
			float a1 = (-2.0f * cs) / a0;
			float a2 = (1.0f - alpha) / a0;
			c1 = a1 + 2.0f;
			c2 = (1.0f + a1 + a2) / c1;
			float bp_b0 = (sn * 0.5f) / a0;
			float bp_d0 = bp_b0;
			float bp_d1 = (2.0f * bp_b0) / c1;
			float bp_d2 = 0.0f;
			float t = lp2bp2hp * 2.0f;
			if (t < 1.0f)// LP -> BP
			{
				float lp_d0 = ((1.0f - cs) * 0.5f) / a0;
				float lp_d1 = c2;
				float lp_d2 = 1.0f;
				d0 = lp_d0 + t * (bp_d0 - lp_d0);
				d1 = lp_d1 + t * (bp_d1 - lp_d1);
				d2 = lp_d2 + t * (bp_d2 - lp_d2);
			}
			else// BP -> HP
			{
				t -= 1.0f;
				float hp_d0 = ((1.0f + cs) * 0.5f) / a0;
				float hp_d1 = 0.0f;
				float hp_d2 = 0.0f;
				d0 = bp_d0 + t * (hp_d0 - bp_d0);
				d1 = bp_d1 + t * (hp_d1 - bp_d1);
				d2 = bp_d2 + t * (hp_d2 - bp_d2);
			}
		}
		void SetFilterParams(float cutoff, float reso, float morph) override
		{
			DesignBasicFilter(cutoff, reso, morph);
		}
		void SetSampleRate(float sr) override
		{
			sampleRate = sr;
		}
	};

	class SVFKernel2order
	{
		float z1 = 0, z2 = 0;
		float d0 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;
	public:
		inline float ProcessSample(float in)
		{
			float x = in - z1 - z2;
			float out = d0 * x + d1 * z1 + d2 * z2;
			z2 += c2 * z1;
			z1 += c1 * x;
			return out;
		}
		void Reset() { z1 = z2 = 0; }
		void SetCoeffs(float* coeffs)//d0,d1,d2,c1,c2
		{
			d0 = coeffs[0];
			d1 = coeffs[1];
			d2 = coeffs[2];
			c1 = coeffs[3];
			c2 = coeffs[4];
		}
	};
	class Elliptic6order :public Filter
	{
	private:
		constexpr static int numSOS = 6 / 2;
		SVFKernel2order svfs[numSOS];
		SVFilter12dB bw12;

		float ellipCutoff = 10000.0;
		//c1,c2, lpd0,lpd1,lpd2, bpd..., hpd...
		constexpr static float ellpCoeffs[3][11] = {
		{ 1.18574364f, 0.358311099f, 0.161611388f, 0.545181548f, 1.52153129f, 0.210616062f, 0.355247213f, 0.0f, 0.274480195f, 0.0f, 0.0f },
		{ 1.43285513f, 0.682587663f, 0.161611388f, 0.451159045f, 0.660953999f, 0.210616062f, 0.293980958f, 0.0f, 0.274480195f, 0.0f, 0.0f },
		{ 1.60344117f, 0.910665734f, 0.161611388f, 0.40316138f, 0.442710607f, 0.210616062f, 0.26270507f, 0.0f, 0.274480195f, 0.0f, 0.0f }
		};

		float sampleRate = 48000.0;
		inline void SVFApplyAPF(float* dst, const float* src, float k)
		{
			const float t0 = k - 1.0;
			const float t1 = t0 * t0;
			const float t2 = k * k * src[3] * src[4];
			const float t3 = k * src[3];
			const float t4 = 1.0 / (-t0 * t3 + t1 + t2);
			const float t5 = src[1] * t0;
			const float t6 = 2 * k * src[4];
			const float t7 = -k + t6 + 1;
			const float t8 = 1.0 / t7;
			const float t9 = k + 1;
			dst[0] = t4 * (src[0] * t1 + src[2] * t2 - t3 * t5);
			dst[1] = t8 * (src[2] * t6 - t5);
			dst[2] = src[2];
			dst[3] = src[3] * t4 * t7 * t9;
			dst[4] = src[4] * t8 * t9;
		}
		void EllipApplyMorph(float* dst, const float* src, float morph)
		{
			const float c1 = src[0];
			const float c2 = src[1];
			const float* lp = src + 2;
			const float* bp = src + 5;
			const float* hp = src + 8;
			float d0, d1, d2;
			float t = morph * 2.0f;
			if (t < 1.0f)
			{
				d0 = lp[0] + t * (bp[0] - lp[0]);
				d1 = lp[1] + t * (bp[1] - lp[1]);
				d2 = lp[2] + t * (bp[2] - lp[2]);
			}
			else
			{
				t -= 1.0f;
				d0 = bp[0] + t * (hp[0] - bp[0]);
				d1 = bp[1] + t * (hp[1] - bp[1]);
				d2 = bp[2] + t * (hp[2] - bp[2]);
			}
			dst[0] = d0;
			dst[1] = d1;
			dst[2] = d2;
			dst[3] = c1;
			dst[4] = c2;
		}
	public:
		void SetFilterParams(float cutoff, float reso, float morph)override
		{
			bw12.SetFilterParams(cutoff, reso, morph);

			if (cutoff > sampleRate / 2 - 3000.0)cutoff = sampleRate / 2 - 3000.0;
			float k =
				cheapSinPi((cutoff - ellipCutoff) / sampleRate) /
				cheapSinPi((cutoff + ellipCutoff) / sampleRate);

			float baseCoeffs[5];
			float warpedCoeffs[5];
			for (int i = 0; i < numSOS; ++i)
			{
				EllipApplyMorph(baseCoeffs, ellpCoeffs[i], morph);
				SVFApplyAPF(warpedCoeffs, baseCoeffs, k);
				svfs[i].SetCoeffs(warpedCoeffs);
			}
		}
		void Reset() override
		{
			bw12.Reset();
			for (int i = 0; i < numSOS; ++i)
				svfs[i].Reset();
		}
		float ProcessSample(float x) override
		{
			x = bw12.ProcessSample(x);
			for (int i = 0; i < numSOS; ++i)
				x = svfs[i].ProcessSample(x);
			return x;
		}
		void SetSampleRate(float sr) override
		{
			bw12.SetSampleRate(sr);
			sampleRate = sr;
		}
	};

	class CombFilter :public Filter
	{
	private:
		constexpr static int MaxBufferSize = 2048;//minfreq = 48000/2048
		float sampleRate = 48000;
		float buf[MaxBufferSize] = { 0 };
		int pos = 0;
		float realDelayTime = 1.0;
		float delayTime = 1.0;
		float fb = 0.0;
		float morph = 0;

		float z = 0;
		inline float ProcessAPF(float x, float k)//k 0->1 : z^-1 -> 1
		{
			x = x - k * z;
			float y = k * x + z;
			z = x;
			return y;
		}
		float ddcz = 0;
		inline float ProcessDeDC(float x)
		{
			ddcz += 0.01 * (x - ddcz);
			return x - ddcz;
		}
	public:
		inline float ProcessSample(float x) override final
		{
			realDelayTime += 0.05 * (delayTime - realDelayTime);
			int idt = realDelayTime;
			float kfrac = realDelayTime - idt;
			int writePos = pos % MaxBufferSize;
			int readPos = (pos + MaxBufferSize - idt) % MaxBufferSize;
			pos++;

			float xWithFB_Add = x + buf[readPos] * fb;
			float xWithFB_Sub = x - buf[readPos] * fb;
			float out_Add = xWithFB_Add + buf[readPos];
			float out_Sub = xWithFB_Sub - buf[readPos];

			float xWithFB = morph < 0.5 ? xWithFB_Add : xWithFB_Sub;
			float out = out_Add + (out_Sub - out_Add) * morph;

			buf[writePos] = ProcessAPF(ProcessDeDC(xWithFB), (1.0 - kfrac) / (1.0 + kfrac));
			return out * 0.5;
		}
		void Reset() override
		{
			for (auto& v : buf)v = 0;
			z = 0;
			ddcz = 0.0f;
			realDelayTime = delayTime;
		}
		void SetSampleRate(float sr) override { sampleRate = sr; }
		void SetFilterParams(float cutoff, float reso, float morph) override
		{
			delayTime = sampleRate / cutoff - 1;
			if (delayTime < 1)delayTime = 1;
			if (delayTime > MaxBufferSize - 1)delayTime = MaxBufferSize - 1;
			fb = 1.0 - 1.0 / reso;
			if (fb < 0)fb = 0;
			if (fb > 0.999)fb = 0.999;
			this->morph = morph;
		}
	};

	class CombFilter4Stage :public Filter
	{
	private:
		CombFilter cf1, cf2, cf3, cf4;
		float gainFix = 1.0;
	public:
		float ProcessSample(float x) override
		{
			x *= gainFix;
			x = cf1.ProcessSample(x);
			x = cf2.ProcessSample(x);
			x = cf3.ProcessSample(x);
			x = cf4.ProcessSample(x);
			return x;
		}
		void Reset() override
		{
			cf1.Reset();
			cf2.Reset();
			cf3.Reset();
			cf4.Reset();
		}
		void SetSampleRate(float sr) override
		{
			cf1.SetSampleRate(sr);
			cf2.SetSampleRate(sr);
			cf3.SetSampleRate(sr);
			cf4.SetSampleRate(sr);
		}
		void SetFilterParams(float cutoff, float reso, float morph) override
		{
			reso = sqrtf(sqrtf(reso));
			float morphGain = 1.0f + 3.0f * cheapSinPi(morph);
			gainFix = morphGain / (reso * reso);
			morph = morph * morph + 1.0;
			cf1.SetFilterParams(cutoff * (morph *= morph), reso, 0);
			cf2.SetFilterParams(cutoff * (morph *= morph), reso, 0);
			cf3.SetFilterParams(cutoff * (morph *= morph), reso, 0);
			cf4.SetFilterParams(cutoff * (morph *= morph), reso, 0);
		}
	};
	class PhaserFilter :public Filter
	{
	private:
		constexpr static int MaxStages = 12;
		int numStages = 2;
		float sampleRate = 48000;

		float z[MaxStages] = { 0 };
		float k = 0, fb = 0, apfout = 0;
		float gainfix = 1.0;
		inline float ProcessAPF(float x)//k 0->1 : z^-1 -> 1
		{
			for (int i = 0; i < numStages; ++i)
			{
				x = x - k * z[i];
				float y = k * x + z[i];
				z[i] = x;
				x = y;
			}
			return x;
		}
		float ddcz = 0;
		inline float ProcessDeDC(float x)
		{
			ddcz += 0.01 * (x - ddcz);
			return x - ddcz;
		}
	public:
		float ProcessSample(float x) override
		{
			float apfin = x * gainfix + fb * apfout;
			apfout = ProcessAPF(ProcessDeDC(apfin));
			return (apfin + apfout) * 0.5;
		}
		void Reset() override
		{
			for (auto& v : z)v = 0;
			apfout = 0;
		}
		void SetSampleRate(float sr) override
		{
			sampleRate = sr;
		}
		void SetFilterParams(float cutoff, float reso, float morph) override
		{
			numStages = morph * (MaxStages - 2) + 2.0;
			float normf = cutoff / sampleRate;
			k = -cheapSinPi(0.5 / numStages - normf) /
				cheapSinPi(0.5 / numStages + normf);
			fb = 1.0 - 1.0 / reso;
			gainfix = sqrtf(fabs(fb - 1.0));
		}
	};
}