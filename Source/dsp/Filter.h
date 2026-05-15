#pragma once
#include <math.h>

namespace MinusMKI
{
	struct SVFilter
	{
		float z1 = 0, z2 = 0;
		float d0 = 0, d1 = 0, d2 = 0, c1 = 0, c2 = 0;
		inline float ProcessSample(float in)
		{
			float x = in - z1 - z2;
			float out = d0 * x + d1 * z1 + d2 * z2;
			z2 += c2 * z1;
			z1 += c1 * x;
			z2 = Nonlinear(z2);
			return out;
		}

		float a = 0.001, I = 1.0;
		const float k = 2.5;
		const float kcoef = (k - logf(1.0 + k)) / (k * k);
		inline float Nonlinear(float x)
		{
			float sgnx = x > 0 ? 1.0 : -1.0;
			x = fabsf(x);
			float f1 = x * a;
			float f2 = x / a;
			float f0 = f1 + (f2 - f1) * (1.0 / (k * x + 1.0));
			f0 = f0 * sgnx / (2.0 * I);
			return f0;
		}
		void DesignNonlinear(float drive)
		{
			a = 1.0 - drive + 0.001;
			I = a * 0.5 + (1.0 / a - a) * kcoef;
		}
		
		void Reset() { z1 = z2 = 0; }
		void DesignBasicFilter(float cutoff, float reso, float lp2bp2hp)
		{
			constexpr float pi = 3.14159265358979323846f;
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