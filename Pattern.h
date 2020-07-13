#ifndef PATTERN_H
#define PATTERN_H

#include "Matrix.h"
#include "Color.h"

namespace rt {


	enum Pattern_type {
		STRIPE,
		GRADIENT,
		RING,
		CHECKER
	};


	class Pattern {
	public:

		Pattern();

		// Pattern type.
		Pattern_type pt;

		// Pattern transform.
		Matrix b_pattern_tx;

		// Inverse Pattern transform;
		Matrix b_inv_pattern_tx;

		// A stripe pattern is a data structure storing a collection of colors
		// that can is defined by the material.
		std::vector<Color> m_stripe_pattern;

		// A pattern is defined.
		bool m_pattern_exist;

		// Computes the appropriate color of a point according to a stripe pattern.
		Color stripe_at(Point& p);

		// Transforms the stripped point.
		Color stripe_at_object( Point& world_p, Matrix& object_itx);

	};

}

#endif // !PATTERN_H
