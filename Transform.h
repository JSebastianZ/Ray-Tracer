/*
* Interface file for the Ray Tracer transform operations.
*
*
* This file contains the declaration for the Tranform class,
* which are used to change points and vectors within a 3D space.
*
*/

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Precision.h"
#include "Point.h"
#include "Vector.h"
#include "Matrix.h"

namespace rt {

	/*
	* Transforms are used to translate, scale, reflect, mirror,
	* rotate and sheare points and vectors.
	*/

	class Transform {

	public:

		// Moves point B according to point P xyz values.
		static Point translate_point(const Point& p, const Point& b);

		// Moves vector V according to point P xyz values.
		static Vector translate_vector(const Point& p, const Vector& v);

		// Scales point B according to point P xyz values.
		static Point scale_point(const Point& p, const Point& b);

		// Scales vector V according to point P xyz values.
		static Vector scale_vector(const Point& p, const Vector& v);

		// Reflects point P in the X axis.
		static Point reflect_point_x(const Point& p);

		// Reflects point P in the Y axis.
		static Point reflect_point_y(const Point& p);

		// Reflects point P in the Z axis.
		static Point reflect_point_z(const Point& p);

		// Reflects point P in the XYZ axis.
		static Point mirror_point(const Point& p);

		// Reflects vector V in the X axis.
		static Vector reflect_vector_x(const Vector& v);

		// Reflects vector V in the Y axis.
		static Vector reflect_vector_y(const Vector& v);

		// Reflects vector V in the Z axis.
		static Vector reflect_vector_z(const Vector& v);

		// Reflects vector V in the XYZ axis.
		static Vector mirror_vector(const Vector& v);

		// Rotates point P around the X axis by D degrees.
		static Point rotate_point_x(const Point& p, real d);

		// Rotates point P around the Y axis by D degrees.
		static Point rotate_point_y(const Point& p, real d);

		// Rotates point P around the Z axis by D degrees.
		static Point rotate_point_z(const Point& p, real d);

		// Rotates vector V around the X axis by D degrees.
		static Vector rotate_vector_x(const Vector& v, real d);

		// Rotates vector V around the Y axis by D degrees.
		static Vector rotate_vector_y(const Vector& v, real d);

		// Rotates vector V around the Z axis by D degrees.
		static Vector rotate_vector_z(const Vector& v, real d);

		// Sheares point P according to the passed parameters.
		static Point shearing(const Point& p, const real& xy, const real& xz, const real& yx,
			const real& yz, const real& zx, const real& zy);
	};
}

#endif // !TRANSFORM_H