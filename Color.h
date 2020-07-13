/*
* Interface file for the Ray Tracer color class.
*
*
* This file contains the definitions for the Color class,
* which provides functions to manipulate RGB color values.
*
*/

#ifndef COLOR_H
#define COLOR_H

#include "Precision.h"
#include "Ray.h"
#include <vector>

namespace rt {

	/*
	* The Ray Tracer camera renders a 3D scene onto a 2D image.
	* It is located one unit away from the final 2D rendered image.
	*/

	class Color{

	public:

		// A color is represented by its RGB values.
		real m_r;
		real m_g;
		real m_b;

		Color();

		// Constructor with RGB values specified.
		Color(real r, real g, real b);

		// Color addition.
		static Color addColors(const Color& x, const Color& y);

		// Color substraction.
		static Color subColors(const Color& x, const Color& y);

		// Modify RGB value by specified number of times to adjust brightness.
		static Color scalarMultiplication(const Color& a, const real& x);

		// Mix two colors.
		static Color multColors(const Color& x, const Color& y);

	};
}

#endif // !COLOR_H
