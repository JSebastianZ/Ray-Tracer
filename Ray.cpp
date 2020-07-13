#include "Ray.h"

using namespace rt;

// Default constructor of a ray at the origin.
Ray::Ray() : m_origin{ Point(0,0,0) }, m_direction{ Vector(0,0,0) }{}

// Construct ray at specified origin and position.
// Ray direction must be normalized.
Ray::Ray(Point origin, Vector position) : m_origin{ origin }, m_direction{ position }{}

// Destructor.
Ray::~Ray() {}

Ray Ray::transform(const Ray& ray, const Matrix& transmatrix) {
	Ray r;
	r.m_origin = Matrix::multiply_point(transmatrix, ray.m_origin);
	// why normalized???
	r.m_direction = Matrix::multiply_vector(transmatrix, ray.m_direction);
	//r.m_direction = Vector::normalize(Matrix::multiply_vector(transmatrix, ray.m_direction));	
	r.m_direction.m_w = 0;
	return r;
	}

// Given a ray, and a hit value, returns the 
// 3D coordinates of the hit point.
Point Ray::position(const Ray& r, const real& d) {
	Vector v(Vector::scalarMultiplication(r.m_direction, d));
	Point p(Point::addVector2Point(r.m_origin, v));
	return p;
}