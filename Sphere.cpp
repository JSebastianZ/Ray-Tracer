#include "Sphere.h"
#include <iostream>
using namespace rt;

// Spheres are defined by their diameter and position within the 3D space.

// Spheres are constructed according to specified parameters.

Sphere::Sphere() : m_diameter{ 2 }, m_center{ Point(0,0,0) }, Shape{ Material(),Matrix::identity(4) }{
	m_radius = (m_diameter / 2);
}

// As objects, a material and transform matrix must be provided upon construction.
Sphere::Sphere(real diameter = 2, Point center = Point(0,0,0), Material material = Material(), Matrix transform = Matrix::identity(4)) : m_diameter{ diameter }, m_center{ center }, Shape{ material,transform } {
	m_radius = (m_diameter / 2);
}

void Sphere::intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap) {
	real t0 = 0;
    real t1 = 0; //intersection points
	real a = Vector::dotProduct(r.m_direction, r.m_direction);
	real b = 2 * (Vector::dotProduct(r.m_direction, Point::subPoints(r.m_origin, this->m_center)));
	real c = (Vector::dotProduct(Point::subPoints(r.m_origin, this->m_center),
		Point::subPoints(r.m_origin, this->m_center)) - (this->m_radius * this->m_radius));
	real discriminant = (pow(b, 2) - (4.0 * a * c));
	if (discriminant > 0) {
		t0 = ((-1 * b) - sqrt(discriminant)) / (2 * a);
		t1 = ((-1 * b) + sqrt(discriminant)) / (2 * a);
		hits.push_back(t0);
		hits.push_back(t1);
		hmap.insert({ t0, this });
		hmap.insert({ t1, this });
	}
	else if ( discriminant == 0) {
		t0 = (-1 * b) / (2 * a);
		t1 = (b) / (2 * a);
		hits.push_back(t0);
		hits.push_back(t1);
		hmap.insert({ t0, this });
		hmap.insert({ t1, this });

	}
	else { hits.push_back(infinite); } // push -999999 

}

void Sphere::normal_at(Point& hitpoint, Matrix& inv_tx, Vector& normal) {
	Point object_point(Matrix::multiply_point(inv_tx, hitpoint));	// Computes normal vector at 3D hit point.
	Vector object_normal(Point::subPoints(object_point, Point(0, 0, 0)));
	Vector world_normal(Matrix::multiply_vector(Matrix::transpose(inv_tx), object_normal));
	world_normal.m_w = 0;
	normal = Vector::normalize(world_normal);
}