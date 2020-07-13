#include "Pattern.h"
#include <iostream>

using namespace rt;

Pattern::Pattern() : m_pattern_exist { true }{};

// Computes the appropriate color of a point according to a stripe pattern.
Color Pattern::stripe_at(Point& p) {

	if (this->pt == STRIPE) {
		int v = floor(p.m_x);
		if ((v % 2) == 0) {
			return this->m_stripe_pattern[0];
		}
		else {
			return this->m_stripe_pattern[1];
		}
	}
	if (this->pt == GRADIENT) {
		Color distance = Color::subColors(this->m_stripe_pattern[0], this->m_stripe_pattern[1]);
		real fraction = p.m_x - floor(p.m_x);
		return Color::addColors(this->m_stripe_pattern[1], Color::scalarMultiplication(distance, fraction));
	}
	if (this->pt == RING) {
		int v = floor(sqrt(pow(p.m_x, 2) + pow(p.m_z, 2)));
		if ((v % 2) == 0) {
			return this->m_stripe_pattern[0];
		}
		else {
			return this->m_stripe_pattern[1];
		}
	}
	if (this->pt == CHECKER) {
		int v = floor(p.m_x) + floor(p.m_y) + floor(p.m_z);
		if ((v % 2) == 0) {
			return this->m_stripe_pattern[0];
		}
		else {
			return this->m_stripe_pattern[1];
		}
	}
}

Color Pattern::stripe_at_object( Point& world_p, Matrix& object_itx) {
	Point object_point = Matrix::multiply_point(object_itx, world_p);
	Point pattern_point = Matrix::multiply_point(this->b_inv_pattern_tx, object_point);
	Color color = stripe_at( pattern_point);
	return color;
}