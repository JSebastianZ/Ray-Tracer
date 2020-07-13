/*
* Interface file for the Ray Tracer objects.
*
*
* This file contains the definitions for the object class,
* which is the base template for all objects within a 3D scene.
*
*/

#ifndef SHAPE_H
#define SHAPE_H

#include "Material.h"
#include "Matrix.h"
#include <stdlib.h>
#include "Ray.h"
#include <vector>
#include <unordered_map>

namespace rt {

	class Shape{
	public:
		// Objects intersected by rays are identified by an ID.
		int b_id;

		// Every shape has a material with color
		// and ambient, diffuse, specular propeties.
		Material b_material;

		// All shapes have a transform matrix
		// that defines their position.
		Matrix b_transform;

		// The inverse matrix is pre-computed to reduce
		// frequent runtime calculations of transformations.
		Matrix b_inv_tx;


		// Default constructor. Assigns a random ID, and an identity tranform matrix.
		Shape();

		// Object constructed with specified material and transform matrix. ID is random.
		Shape(Material m, Matrix t);

		// Destructor.
		~Shape();

		// Returns an Object ID.
		int getID() const;
	public:
		// Ray intersection function.
		virtual void intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap) = 0;

		// Normal at hipoint.
		virtual void normal_at(Point& hitpoint, Matrix& inv_tx , Vector& normal) = 0;


	};
}

#endif // !SHAPE_H
