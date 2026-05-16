#pragma once
#include <math.h>

namespace MinusMKI
{
	class Filter
	{

	};
	class SVFilter
	{
	private:
		float z1 = 0, z2 = 0;
		float d0 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;

		float a = 0.25, d = 1.0;
	public:
		inline float ProcessSample(float in)
		{
			float x = in - z1 - z2;
			float out = d0 * x + d1 * z1 + d2 * z2;
			z2 += c2 * z1;
			z1 += c1 * x;
			return Nonlinear(out);
		}
		inline float Nonlinear(float x)
		{
			bool sgnx = x >= 0;
			x = fabsf(x);
			float y = x + (a - x) * (1.0 - 1.0 / (d * x + 1.0));
			return sgnx ? y : -y;
		}
		void DesignNonlinear(float drive)
		{
			d = drive * drive * 100.0;
		}

		void Reset() { z1 = z2 = 0; }
		void DesignBasicFilter(float cutoff, float reso, float lp2bp2hp, float sampleRate = 48000.0)
		{
			constexpr float pi = 3.14159265358979323846f;

			if (cutoff < 40.0)reso = 0.707;
			cutoff /= sampleRate;
			if (cutoff < 0.00005)cutoff = 0.00005;
			if (cutoff > 0.45)cutoff = 0.45;

			float w0 = 2.0f * pi * cutoff;
			float cs = cosf(w0);
			float sn = sinf(w0);
			float alpha = sn / (2.0f * reso); // reso = Q
			float a0 = 1.0f + alpha;
			float b0 = ((1.0f - cs) * 0.5f) / a0;
			//float b1 = (1.0f - cs) / a0;
			//float b2 = ((1.0f - cs) * 0.5f) / a0;
			float a1 = (-2.0f * cs) / a0;
			float a2 = (1.0f - alpha) / a0;

			c1 = a1 + 2.0f;
			c2 = (1.0f + a1 + a2) / c1;
			// LPF cookbook -> SVF
			float lp_d0 = ((1.0f - cs) * 0.5f) / a0;
			float lp_d1 = c2;
			float lp_d2 = 1.0f;
			// BPF cookbook, constant skirt gain, peak gain = Q
			float bp_b0 = (sn * 0.5f) / a0;
			float bp_d0 = bp_b0;
			float bp_d1 = (2.0f * bp_b0) / c1;
			float bp_d2 = 0.0f;
			// HPF cookbook -> SVF
			float hp_d0 = ((1.0f + cs) * 0.5f) / a0;
			float hp_d1 = 0.0f;
			float hp_d2 = 0.0f;
			float t = lp2bp2hp * 2.0f;
			if (t < 1.0f)// LP -> BP
			{
				d0 = lp_d0 + t * (bp_d0 - lp_d0);
				d1 = lp_d1 + t * (bp_d1 - lp_d1);
				d2 = lp_d2 + t * (bp_d2 - lp_d2);
			}
			else// BP -> HP
			{
				t -= 1.0f;
				d0 = bp_d0 + t * (hp_d0 - bp_d0);
				d1 = bp_d1 + t * (hp_d1 - bp_d1);
				d2 = bp_d2 + t * (hp_d2 - bp_d2);
			}
		}
	};
}