/*
* Interface file for the Ray Tracer tuples.
*
*
* This file contains the definitions for the Tuple class,
* which is the base class for Points and Vectors.
*
*/

#ifndef TUPLE_H
#define TUPLE_H

#include "Precision.h"

namespace rt {

	/*
	* Tuples are used to represent  3D primitives in Cartesian plane (e.g. point, vectors).
	*/

	class Tuple {

	public:

		// Holds the value along the x axis.
		real m_x;

		// Holds the value along the y axis.
		real m_y;

		// Holds the value along the z axis.
		real m_z;

		// Indicates primitive type. 0 = vector / 1 = point.
		int m_w;

		// Default constructor.
		Tuple();

		// Constructs a tuple according to parameters passed.
		Tuple(real x, real y, real z, int w);

		// Destructor.
		~Tuple();

	};
}

#endif // !TUPLE_H
