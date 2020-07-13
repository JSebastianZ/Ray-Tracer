#include "Camera.h"

using namespace rt;

// Default constructor.
Camera::Camera() : m_hsize{ 0 }, m_vsize{ 0 }, m_field_of_view{ 0 }, m_transform{ Matrix::identity(4) } {}

// Constructs camera with specified parameter values, and identity matrix as transform.
Camera::Camera(real hsize, real vsize, real field_of_view) : m_hsize{ hsize }, m_vsize{ vsize }, m_field_of_view{ field_of_view }, m_transform{ Matrix::identity(4) }{ m_inv_tx = Matrix::inverse(m_transform); }

// Constructor for camera object.
Camera::Camera(real hsize, real vsize, real field_of_view, Matrix transform) : m_hsize{ hsize }, m_vsize{ vsize }, m_field_of_view{ field_of_view }, m_transform{ transform }{m_inv_tx = Matrix::inverse(m_transform); }


real Camera::pixel_size() {
	this->half_view = std::tan(this->m_field_of_view / 2.0);
	this->aspect = (this->m_hsize / this->m_vsize);

	if (aspect >= 1) {
		this->m_half_width = this->half_view;
		this->m_half_height = (this->half_view / this->aspect);
	}
	else {
		this->m_half_width = (this->half_view * this->aspect);
		this->m_half_height = this->half_view;
	}
	return ((this->m_half_width * 2.0) / this->m_hsize);
}

Ray Camera::ray_for_pixel(const int& x, const int& y) {
	real xoffset = (x + 0.5) * this->pixel_size();			// Offset from edge of canvas to pixels center
	real yoffset = (y + 0.5) * this->pixel_size();
	real world_x = this->m_half_width - xoffset;			// Untransformed coordinates of pixel in world space
	real world_y = this->m_half_height - yoffset;
	Point current = (Point(world_x, world_y, -1));
	Point t_origin = (Point(0, 0, 0));
	Point pixel = Matrix::multiply_point(this->m_inv_tx, current);
	Point origin = Matrix::multiply_point(this->m_inv_tx, t_origin);
	Vector direction = Vector::normalize(Point::subPoints(pixel, origin));
	Ray ray = Ray(origin, direction);
	return ray;
}
