#include "Color.h"

using namespace rt;

real m_r = 0;
real m_g = 0;
real m_b = 0;

Color::Color() : m_r{ 0 }, m_g{ 0 }, m_b{ 0 }{}

// Constructor with RGB values specified.
Color::Color(real r, real g, real b) : m_r{ r }, m_g{ g }, m_b{ b } {}

// Color addition.
Color Color::addColors(const Color& x, const Color& y) {
	real r = x.m_r + y.m_r;
	real g = x.m_g + y.m_g;
	real b = x.m_b + y.m_b;
	return{ r, g, b };
}

// Color substraction.
Color Color::subColors(const Color& x, const Color& y) {
	real r = x.m_r - y.m_r;
	real g = x.m_g - y.m_g;
	real b = x.m_b - y.m_b;
	return{ r, g, b };
}

// Modify RGB value by specified number of times to adjust brightness.
Color Color::scalarMultiplication(const Color& a, const real& x) {
	real r = a.m_r * x;
	real g = a.m_g * x;
	real b = a.m_b * x;
	return{ r, g, b };
}

// Mix two colors.
Color Color::multColors(const Color& x, const Color& y) {
	real r = x.m_r * y.m_r;
	real g = x.m_g * y.m_g;
	real b = x.m_b * y.m_b;
	return{ r, g, b };
}

