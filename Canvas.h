/*
* Interface file for the Ray Tracer 2D rendered image canvas.
*
*
* This file contains the definitions for the canvas class,
* which provides the a plane where images will be draw.
*
*/

#ifndef CANVAS_H
#define CANVAS_H

#include "Color.h"

namespace rt {
	class Canvas {

		public:

		// Canvas width and heigth;
		int m_height;
		int m_width;

		// A canvas is implemented as a one dimention array of pixels, each one storing an RGB color value.
		// Canvas is composed of colors for each pixels, not composition by aggregation of Colors!
		Color* pixels;

		// Canvas constructor. The pixels array size is width * height
		Canvas(const int width, const int height);

		// Draws RGB color values into a single canvas pixel.
		void write_pixel(int x, int y, Color color);

		// Returns the specified canvas pixel RGB color value.
		Color get_pixel(int x, int y);

	};
}

#endif // !CANVAS_H
