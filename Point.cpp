#include "Point.h"

using namespace rt;

// Default constructor. Point in the origin of 3D space.
Point::Point() : Tuple{ 0,0,0,1 } {}

// Creates a point object at the specified xyz position.
Point::Point(real x, real y, real z) : Tuple{x,y,z,1} {}

// Destructor.
Point::~Point() {}

// Returns the Point resulting from the addition of points A and B.
Point Point::addPoints(const Point& a, const Point& b) {
	real x = a.m_x + b.m_x;
	real y = a.m_y + b.m_y;
	real z = a.m_z + b.m_z;
	return { x, y, z };
}

// Returns the Point resulting from the addition of point A and vector B.
Point Point::addVector2Point(const Point& a, const Tuple& b) {
	real x = a.m_x + b.m_x;
	real y = a.m_y + b.m_y;
	real z = a.m_z + b.m_z;
	return { x, y, z };
}

// Returns the Vector resulting from substracting point B from A.
// A - B returns a vector point from B towards A.
Tuple Point::subPoints(const Point& a, const Point& b) {
	real x = a.m_x - b.m_x;
	real y = a.m_y - b.m_y;
	real z = a.m_z - b.m_z;
	int w = 0;
	return { x, y, z, w };
}

// Returns the point resulting from substracting vector B from point A.
Point Point::subVectorFromPoint(const Point& a, const Tuple& b) {
	real x = a.m_x - b.m_x;
	real y = a.m_y - b.m_y;
	real z = a.m_z - b.m_z;
	return { x, y, z };
}

// Returns the point that results from scaling point A by I.
Point Point::scalarMultiplication(const Point& a, const real& i) {
	return { i * a.m_x, i * a.m_y, i * a.m_z };
}

// Returns the point that results from dividing point A by I.
Point Point::scalarDivision(const Point& a, const real& i) {
	return { a.m_x / i, a.m_y / i, a.m_z / i };
}