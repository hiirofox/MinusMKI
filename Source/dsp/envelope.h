#pragma once
#include <math.h>

namespace MinusMKI
{
	enum class EnvelopeMode {
		PolyphonicResetOnNoteOn = 0,	//复音，按下时重置
		PolyphonicNoReset = 1,			//复音，无重置
		GlobalNoReset = 2,				//全局，无重置
		GlobalResetOnFirstNoteOn = 3,	//全局，按下第一个键重置
		GlobalResetOnAnyNoteOn = 4		//全局，按下任意键重置
	};
	class Envelope
	{
	private:
	public:
		virtual void Reset();
		virtual void SetSampleRate(float sr) {};
		virtual void SetNoteState(bool off0on1) {};
		virtual void SetParams(float p1, float p2, float p3, float p4, float p5, float p6) {};
		virtual void Step() {};
		virtual float GetValue() {};
	};
	class ADSR :public Envelope
	{
	private:
		float sampleRate = 48000;
		//0:release 1:attack 2:decay 3:sustain
		float s[4] = { 0 }, k[4] = { 0 }, r[4] = { 0 };
		float pos = 0;
		float y = 0;
	public:
		void Step() override
		{
			int idx = pos;
			y = s[idx] * y + k[idx];
			pos += r[idx];
		}
		void Reset() override
		{
			pos = 0.0;
			y = 0.0;
		}
		void SetNoteState(bool state) override
		{
			if (state)
			{
				pos = 1.0;
				y = 0.0;
			}
			else
			{
				pos = 0.0;
			}
		}
		float GetValue() override { return y; }
		void SetSampleRate(float sr) override { sampleRate = sr; }
		void SetSegment(float y0, float y1, float len, float shape, int i)
		{
			if (len < 4.0)len = 4.0;
			float sv = 0, kv = 0, rv = 0;
			if (fabsf(shape) > 0.01)
			{
				sv = expf(shape / len);
				float expshape = expf(shape);
				kv = (y1 - y0 * expshape) * (sv - 1.0) / (expshape - 1.0);
				rv = 1.0 / len;
			}
			else
			{
				sv = 1.0;
				kv = (y1 - y0) / len;
				rv = 1.0 / len;
			}
			s[i] = sv;
			k[i] = kv;
			r[i] = rv;
		}
		float EvalSegment(float y0, float y1, float shape, float phase)
		{
			if (fabsf(shape) > 0.01f)
			{
				float expshape = expf(shape);
				float expphase = expf(shape * phase);
				return y0 * expphase + (y1 - y0 * expshape) * (expphase - 1.0f) / (expshape - 1.0f);
			}
			return y0 + (y1 - y0) * phase;
		}
		void SetParams(float attMs, float attShape, float decMs, float decShape, float susV, float relMs) override
		{
			float attl = sampleRate / 1000.0 * attMs;
			float decl = sampleRate / 1000.0 * decMs;
			float rell = sampleRate / 1000.0 * relMs;
			//release
			if (rell < 4.0)rell = 4.0;
			s[0] = expf(-6.90775527898f / rell);//rt60
			k[0] = 0.0;
			r[0] = 0.0;
			//sustain
			s[3] = 1.0;
			k[3] = 0.0;
			r[3] = 0.0;
			//attack
			SetSegment(0.0, 1.0, attl, attShape, 1);
			//decay
			SetSegment(1.0, susV, decl, decShape, 2);
			//update y
			int idx = pos;
			float posFrac = pos - idx;
			if (idx == 1) y = EvalSegment(0.0f, 1.0f, attShape, posFrac);
			else if (idx == 2) y = EvalSegment(1.0f, susV, decShape, posFrac);
			else if (idx == 3) y = susV;
		}
	};
}