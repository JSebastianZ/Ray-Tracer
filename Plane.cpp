#include "Plane.h"
#include <iostream>
using namespace rt;

// Default plane is over XZ plane, and a point over this plane is (1,0,1).
Plane::Plane() : m_plane_normal{ Vector(0,1,0) }, m_plane_distance{Point(1,0,1)} {}

void Plane::intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap) {
	real c = Vector::dotProduct(r.m_direction, this->m_plane_normal);
	if (abs(c) <= epsilon) {
		hits.push_back(infinite);
	}
	else {
		real t0 = (Vector::dotProduct(Point::subPoints(this->m_plane_distance, r.m_origin), this->m_plane_normal))
			/ (Vector::dotProduct(r.m_direction, this->m_plane_normal));

		hits.push_back(t0);
		hmap.insert({ t0,this });
	}
}

void Plane::normal_at(Point& hitpoint, Matrix& inv_tx, Vector& normal) {
	/*	// Pending formal normal calculator.
	Point object_point(Matrix::multiply_point(inv_tx, hitpoint));	// Computes normal vector at 3D hit point.
	Vector object_normal(Point::subPoints(object_point, Point(0, 0, 0)));
	Vector world_normal(Matrix::multiply_vector(Matrix::transpose(inv_tx), object_normal));
	world_normal.m_w = 0;
	normal = Vector::normalize(world_normal);
	*/
	normal = this->m_plane_normal;
}