/*
* Interface file for the Ray Tracer Sphere object.
*
*
* This file contains the definitions for the Sphere class,
* which provides the interfaces to create sphere objects.
*
*/

#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"
#include <vector>
#include <unordered_map>

namespace rt {

	/*
	* The Ray Tracer spheres derive from the Object class.
	*/



	class Sphere : public Shape {
	public:

		// Spheres are defined by their diameter and position within the 3D space.
		real m_diameter;
		real m_radius;
		Point m_center;

		// Spheres are constructed according to specified parameters.
		// Default constructor creates an sphere centered at the origin,
		// with default material, diameter 2, and identity matrix as transform.
		Sphere();

		// As objects, a material and transform matrix must be provided upon construction.
		Sphere(real diameter, Point center, Material material, Matrix transform);

	// Sphere intersection computation function.
	 void intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap);

	 // Sphere normal computation function at hitpoint.
	 void normal_at(Point& hitpoint, Matrix& inv_tx, Vector& normal);
	};
}

#endif // !SPHERE_H
