/*
* Interface file for the Ray Tracer Plane object.
*
*
* This file contains the definitions for the Sphere class,
* which provides the interfaces to create sphere objects.
*
*/

#ifndef PLANE_H
#define PLANE_H

#include "Shape.h"
#include <vector>
#include <unordered_map>

namespace rt {

	/*
	* The Ray Tracer spheres derive from the Object class.
	*/


	class Plane : public Shape {
	public:

		// A plane has a normal vector perpendicular to its surface.
		Vector m_plane_normal;

		// Point that represents how far is the plane from the world origin.
		Point m_plane_distance;

		// Planes are constructed according to specified parameters.
		// Default constructor creates a plane centered at the origin,
		// with default material, and identity matrix as transform.
		Plane();

		// Plane intersection computation function.
		void intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap);

		// Plane normal computation function at hitpoint.
		void normal_at(Point& hitpoint, Matrix& inv_tx, Vector& normal);

	};

}

#endif // !PLANE_H
