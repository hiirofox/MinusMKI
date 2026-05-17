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
}