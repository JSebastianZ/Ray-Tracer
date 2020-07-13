#ifndef PRECISION_H
#define PRECISION_H

#include <algorithm>
#include <cmath>

namespace rt {
		/*
			Defines a real number precision. RTC can be compiled in
			single or double pricision, just by changing the typedef.
		*/

	typedef double real;

	// Some constants
	constexpr double pi{ 3.14159265 };
	constexpr int infinite{ -999999 };
	constexpr int high_number{ 9999 };

	// Test two for real numbers for equality.
	constexpr double epsilon{ 0.00001 };
	bool approxEqual(double a, double b, double relEpsilon);
}

#endif // !PRECISION_H
