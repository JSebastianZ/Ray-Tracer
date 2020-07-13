/*
* Interface file for the Ray Tracer vectors.
*
*
* This file contains the definitions for the Vector class,
* which are used to represent rays direction, among others.
*
*/

#ifndef VECTOR_H
#define VECTOR_H

#include "Tuple.h"

namespace rt {

	/*
	* Vectors are defined in 3D space by their x, y, and z values.
	* This class defines their construction and operations.
	*/

	class Vector : public Tuple {

	public:

		// 
		// Vectors are implemented as tuples holding x, y, and z values.
		// m_w member is '0' to distinguish from other tuples (e.g. points).

		// Default constructor places a vector in the origin.
		Vector();

		// Create Vector from Tuple
		Vector(Tuple t);

		// Construct vector according to passed parameters.
		Vector(real x, real y, real z);

		// Destructor
		~Vector();

		// Returns a Vector resulting from the addition of two vectors.
		// A + B
		static Vector addVectors(const Vector& a, const Vector& b);

		// Returns a Vector resulting from substracting vector B from A.
		// A - B
		static Vector subVectors(const Vector& a, const Vector& b);

		// Returns the vector result from negating all vector A xyz values.
		static Vector negate(const Vector& a);

		// Returns vector resulting from mutiplying vector A by X.
		static Vector scalarMultiplication(const Vector& a, const real& x);

		// Returns vector resulting from dividing vector A by X.
		Vector scalarDivision(const Vector& a, const real& x);

		// Returns the magnitude of vector A
		// |A|
		static real magnitude(const Vector& a);

		// Returns the vector A normalized.
		static Vector normalize(const Vector& a);

		// Returns the dot product of vectors A and B.
		// A * B
		static real dotProduct(const Vector& a, const Vector& b);

		// Returns the cross product of vectors A and B.
		// A x B
		static Vector crossProduct(const Vector& a, const Vector& b);

		// Returns the reflected vector resulting from the IN vector 
		// intersecting an object at a point with NORMAL vector.
		static Vector reflected(Vector& in, Vector& normal);

	};
}
#endif // !VECTOR_H