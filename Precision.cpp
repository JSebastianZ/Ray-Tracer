#include "Precision.h"

using namespace rt;

	/*
	* This function is implemented to compare two numbers and
	* return TRUE if their difference is smaller than EPSILON.
	*/

	bool approxEqual(double a, double b, double relEpsilon) {
		// Return true if the numbers are close to each other and difference is lower than EPSILON.
		double diff{ std::abs(a - b) };
		if (diff <= epsilon)
			return true;
		// Return true if the difference between a and b is within epsilon % of the larger of a and b (D. Knuth).
		return (diff <= (std::max(std::abs(a), std::abs(b)) * relEpsilon));
	}
