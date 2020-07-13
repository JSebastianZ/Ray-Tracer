/*
* Interface file for the Ray Tracer camera.
*
*
* This file contains the definitions for the camera class,
* which provides the "view" of the rendered scene.
*
*/

#ifndef CAMERA_H
#define CAMERA_H

#include "Precision.h"
#include "Matrix.h"
#include <math.h>
#include "Ray.h"

namespace rt {

	/*
	* The Ray Tracer camera renders a 3D scene onto a 2D image.
	* It is located one unit away from the final 2D rendered image.
	*/

	class Camera {

	public:

		// Horizontal size (in pixels) of the final 2D rendered image.
		real m_hsize;

		// Vertical size (in pixels) of the final 2D rendered image.
		real m_vsize;

		// Angle (radians) of aperture of the camera. Large angle translates into a "zoom out" view.
		real m_field_of_view;

		// Aspect is the ratio of the horizontal size of the final 2D image, to its vertical size.
		real aspect;

		// The width of half of the 2D rendered image. If >= 1, then it is half the horizontal size.
		// If < 1, then it is half the vertical size.
		real half_view;

		// Hold half of the 2D rendered image width and heigth, necessary for pixel size computation.
		real m_half_width;
		real m_half_height;

		// Describes how a "scene" moves relative to the camera.
		Matrix m_transform;

		// Inverse of camera transform, used to reduce computations on runtime.
		Matrix m_inv_tx;

		// Default constructor.
		Camera();

		// Constructs camera with specified parameter values, and identity matrix as transform.
		Camera(real hsize, real vsize, real field_of_view);

		// Constructor for camera object.
		Camera(real hsize, real vsize, real field_of_view, Matrix transform);

		// Destructor
		~Camera() {}

		// Computes the size of a pixel of the final 2D rendered image
		// by dividing its width by its horizontal size.
		// The returned value is used by Ray::ray_for_pixel in order to
		// determine the exact point where a ray must pass through the pixel.
		real pixel_size();

		// Computes 3D coordinates at the center of a given pixel
		// and constructs a ray that passes through it.
		Ray ray_for_pixel(const int& x, const int& y);

	};
}

#endif // !CAMERA_H