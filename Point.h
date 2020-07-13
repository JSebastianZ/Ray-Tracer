/*
* Interface file for the Ray Tracer Point class.
*
*
* This file contains the definitions for the Point class,
* which provides the "view" of the rendered scene.
*
*/

#ifndef POINT_H
#define POINT_H

#include "Tuple.h"

namespace rt {

		/*
		* Points in the Ray Tracer are implemented as tuples of 4 elementes,
		* corresponding to X, Y, and Z position values within the 3D space.
		* The last "w" element value is always "1", and is used to
		* differentiate a point from a vector (tuple with w = 0).
		*/

		class Point : public Tuple {

		public:

			// Points consist of X,Y, and Z values. Points m_w member is always "1". 

			// Default constructor. Point in the origin of 3D space.
			Point();

			// Creates a point object at the specified xyz position.
			Point(real x, real y, real z);

			// Destructor.
			~Point();

			// Returns the Point resulting from the addition of points A and B.
			static Point addPoints(const Point& a, const Point& b);

			// Returns the Point resulting from the addition of point A and vector B.
			static Point addVector2Point(const Point& a, const Tuple& b);

			// Returns the Vector resulting from substracting point B from A.
			// A - B returns a vector point from B towards A.
			static Tuple subPoints(const Point& a, const Point& b);

			// Returns the point resulting from substracting vector B from point A.
			static Point subVectorFromPoint(const Point& a, const Tuple& b);

			// Returns the point that results from scaling point A by I.
			static Point scalarMultiplication(const Point& a, const real& i);

			// Returns the point that results from dividing point A by I.
			static Point scalarDivision(const Point& a, const real& i);
	};
}

#endif // !POINT_H
