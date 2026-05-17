#pragma once
#include <math.h>

namespace MinusMKI
{
	class Filter
	{
	protected:
		inline float cheapCosPi(float x)
		{
			return 1.0 + x * x * (4.0 * x - 6.0);
		}
		inline float cheapSinPi(float x)
		{
			float t = x * (1.0 - x);
			return t * (3.14159265 + 3.40185714 * t);
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

	template<int Order>
	class SVFKernel
	{
	private:
		float z[Order + 1] = { 0 };
		float zeros[Order + 1] = { 0 };
		float* c = zeros;
		float* d = zeros;
	public:
		void Reset()
		{
			for (int i = 0; i <= Order; ++i)z[i] = 0;
		}
		float ProcessSample(float in)
		{
			float x = in;
			for (int i = 1; i <= Order; ++i)
				x -= z[i];
			float out = d[0] * x;
			for (int i = 1; i <= Order; ++i)
				out += d[i] * z[i];
			for (int i = Order; i >= 2; --i)
				z[i] += c[i] * z[i - 1];
			z[1] += c[1] * x;
			return out;
		}
		void SetCDCoeffsData(float* c, float* d)
		{
			this->c = c;
			this->d = d;
		}
	};
	class Elliptic6order :public Filter
	{
	private:
		//0.1dB -60dB 1000Hz@48000Hz
		const float ellipCutoff = 10000.0;
		const double b[7] = { 0.023283907371774663, 0.070294830598407673, 0.12911548236576503, 0.15443551523696455, 0.12911548236576506, 0.070294830598407659, 0.023283907371774659 };
		const double a[7] = { 1, -1.7779600658138057, 2.6506815135932045, -2.2572291693627999, 1.4278113364993401, -0.54822749048802277, 0.11169346554238167 };
		float c[7] = { 0 };
		float d[7] = { 0 };
		SVFKernel<6> svf;
		float sampleRate = 48000;
	public:
		void SetFilterParams(float cutoff, float reso, float morph)override
		{
			if (cutoff > sampleRate / 2 - 3000.0)cutoff = sampleRate / 2 - 3000.0;
			float k = sin(M_PI * (cutoff - ellipCutoff) / sampleRate) / sin(M_PI * (cutoff + ellipCutoff) / sampleRate);

			const float t0 = k * a[1];
			const float t1 = k * k;
			const float t2 = t1 * a[2];
			const float t3 = k * k * k;
			const float t4 = t3 * a[3];
			const float t5 = t1 * t1;
			const float t6 = t5 * a[4];
			const float t7 = t5 * k;
			const float t8 = t7 * a[5];
			const float t9 = t3 * t3;
			const float t10 = 1.0 / (t0 + t2 + t4 + t6 + t8 + t9 * a[6] + 1);
			const float t11 = k * b[1];
			const float t12 = t1 * b[2];
			const float t13 = t3 * b[3];
			const float t14 = t5 * b[4];
			const float t15 = t7 * b[5];
			const float t16 = 6 * t7;
			const float t17 = 2 * k;
			const float t18 = t17 * a[2] + 6;
			const float t19 = 5 * t5;
			const float t20 = t19 * a[5] + 3 * t4;
			const float t21 = 3 * a[3];
			const float t22 = 5 * t0 + t1 * t21;
			const float t23 = 4 * t3;
			const float t24 = 4 * t2 + t23 * a[4];
			const float t25 = t16 * a[6] + t18 + t20 + t22 + t24 + 2 * t6 + t8 + a[1];
			const float t26 = 1.0 / t25;
			const float t27 = t17 * b[2] + 6 * b[0];
			const float t28 = 3 * t13 + t19 * b[5];
			const float t29 = 3 * b[3];
			const float t30 = t1 * t29 + 5 * t11;
			const float t31 = 4 * t12 + t23 * b[4];
			const float t32 = 8 * t3;
			const float t33 = 15 * t5;
			const float t34 = k * t21 + 5 * a[1];
			const float t35 = 6 * t1;
			const float t36 = 8 * k;
			const float t37 = t35 * a[4] + t36 * a[2] + 15;
			const float t38 = 9 * t1;
			const float t39 = 10 * t3;
			const float t40 = 10 * t0 + t38 * a[3] + t39 * a[5];
			const float t41 = 6 * t2 + t20 + t32 * a[4] + t33 * a[6] + t34 + t37 + t40 + t6 + a[2];
			const float t42 = 1.0 / t41;
			const float t43 = k * t29 + 5 * b[1];
			const float t44 = t35 * b[4] + t36 * b[2] + 15 * b[0];
			const float t45 = 10 * t11 + t38 * b[3] + t39 * b[5];
			const float t46 = 12 * k;
			const float t47 = 12 * t1;
			const float t48 = 20 * t3;
			const float t49 = 4 * k;
			const float t50 = t49 * a[4] + 4 * a[2];
			const float t51 = 9 * k;
			const float t52 = 10 * t1;
			const float t53 = t51 * a[3] + t52 * a[5] + 10 * a[1];
			const float t54 = t24 + t4 + t40 + t46 * a[2] + t47 * a[4] + t48 * a[6] + t50 + t53 + a[3] + 20;
			const float t55 = 1.0 / t54;
			const float t56 = t49 * b[4] + 4 * b[2];
			const float t57 = t51 * b[3] + t52 * b[5] + 10 * b[1];
			const float t58 = 15 * t1;
			const float t59 = 5 * k;
			const float t60 = t21 + t59 * a[5];
			const float t61 = t2 + t22 + t36 * a[4] + t37 + t53 + t58 * a[6] + t60 + 6 * a[2] + a[4];
			const float t62 = 1.0 / t61;
			const float t63 = t29 + t59 * b[5];
			const float t64 = 6 * k;
			const float t65 = t0 + t18 + t34 + t50 + t60 + t64 * a[6] + 2 * a[4] + a[5];
			const float t66 = 1.0 / t65;
			const float t67 = a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + 1;
			const float t68 = k + 1;
			c[0] = 0;
			d[0] = t10 * (t11 + t12 + t13 + t14 + t15 + t9 * b[6] + b[0]);
			d[1] = t26 * (2 * t14 + t15 + t16 * b[6] + t27 + t28 + t30 + t31 + b[1]);
			d[2] = t42 * (6 * t12 + t14 + t28 + t32 * b[4] + t33 * b[6] + t43 + t44 + t45 + b[2]);
			d[3] = t55 * (t13 + t31 + t45 + t46 * b[2] + t47 * b[4] + t48 * b[6] + t56 + t57 + 20 * b[0] + b[3]);
			d[4] = t62 * (t12 + t30 + t36 * b[4] + t44 + t57 + t58 * b[6] + t63 + 6 * b[2] + b[4]);
			d[5] = t66 * (t11 + t27 + t43 + t56 + t63 + t64 * b[6] + 2 * b[4] + b[5]);
			d[6] = (b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6]) / t67;
			c[1] = t10 * t25 * t68;
			c[2] = t26 * t41 * t68;
			c[3] = t42 * t54 * t68;
			c[4] = t55 * t61 * t68;
			c[5] = t62 * t65 * t68;
			c[6] = t66 * t67 * t68;
			svf.SetCDCoeffsData(c, d);
		}
		void Reset() override { svf.Reset(); }
		float ProcessSample(float x) override { return svf.ProcessSample(x); }
		void SetSampleRate(float sr) override { sampleRate = sr; }
	};
}