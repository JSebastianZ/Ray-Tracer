#ifndef LIGHT_H
#define LIGHT_H

#include "Point.h"
#include "Color.h"

namespace rt {
	class Light {
	public:

		// Light source position within a 3D scene.
		Point m_source;

		// The intensity of the light source is represented by its RGB color.
		Color m_intensity;

		// Default constructor with light source at the origin and no illumination.
		Light();

		// Light object constructor with specified position and color.
		Light(Point source, Color intensity);

	};
}

#endif // !LIGHT_H
