/*
* Interface file for the Ray Tracer 3D scene.
*
*
* This file contains the definitions for the world class,
* which defines all the objects to be rendered by the Ray Tracer.
*
*/

#ifndef WORLD_H
#define WORLD_H

#include "Sphere.h"
#include "Light.h"
#include "Plane.h"

namespace rt {
	class World {
	public:

		// Every 3D scene must have at least one light source.
		Light m_light_source;

		// World objects are stored in a collection, implemented as a vector.
		std::vector<Shape*> m_objects;

		// Default World sample
		World();

	};
}

#endif // !WORLD_H
