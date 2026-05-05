#pragma once

class DownSampler2x
{
private:
	const float b[5] = { 3.54886431e-04, -7.21089910e-05,  4.74055716e-04 ,-7.21089910e-05, 3.54886431e-04 };
	const float a[5] = { 1.0 ,-3.61096661,  4.93830966, -3.02874037 , 0.70244897 };
	float ProcessFilter(float x)
	{
		const float y = b[0] * x + z0;
		z0 = b[1] * x - a[1] * y + z1;
		z1 = b[2] * x - a[2] * y + z2;
		z2 = b[3] * x - a[3] * y + z3;
		z3 = b[4] * x - a[4] * y;
		return y;
	}
	float out = 0;
	float z0 = 0, z1 = 0, z2 = 0, z3 = 0;
public:
	DownSampler2x()
	{
		Reset();
	}
	void Reset()
	{
		out = 0;
		z0 = 0, z1 = 0, z2 = 0, z3 = 0;
	}
	void ProcessIn(float x)
	{
		out = ProcessFilter(x);
	}
	float GetProcessOut() const
	{
		return out;
	}
};
