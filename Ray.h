/*
* Interface file for the Ray Tracer rays.
*
*
* This file contains the definitions for the ray class,
* which create rays and find intersections with objects.
*
*/

#ifndef RAY_H
#define RAY_H

#include "Point.h"
#include "Vector.h"
#include "Matrix.h"

namespace rt {

	class Ray {
	public:

		// Every ray has an origin and direction
		Point m_origin;
		Vector m_direction;

		// m_inside is used to identify when a ray or eye vector
		// originates inside an object's volume.
		bool m_inside = false;

		// Point in the 3D world where a ray intersects an object.
		Point m_choke;

		// m_hit is the magnitude of m_choke;
		real m_hit;

		// Default constructor of a ray at the origin.
		Ray();

		// Construct ray at specified origin and position.
		// Ray direction must be normalized.
		Ray(Point origin, Vector position);

		// Destructor.
		~Ray();

		// Apply the inverse of the transformation matrix
		// of an Object to a Ray. Returns a new ray with
		// transformed origin and direction. This is used 
		// to render transformed objects.
		static Ray transform(const Ray& r, const Matrix& transmatrix);

		// Given a ray, and a hit value, returns the 
		// 3D coordinates of the hit point.
		static Point position(const Ray& r, const real& d);

	};

}

#endif // !RAY_H
