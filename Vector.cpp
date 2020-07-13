#include "Vector.h"

using namespace rt;

// Default constructor places a vector in the origin.
Vector::Vector() : Tuple{ 0,0,0,0 } {}

// Create Vector from Tuple
Vector::Vector(Tuple t) : Tuple{ t.m_x,t.m_y,t.m_z,0 } {}

// Construct vector according to passed parameters.
Vector::Vector(real x, real y, real z) : Tuple{ x,y,z,0 } {}

// Destructor
Vector::~Vector() {}

// Returns a Vector resulting from the addition of two vectors.
// A + B
Vector Vector::addVectors(const Vector& a, const Vector& b) {
	real x = a.m_x + b.m_x;
	real y = a.m_y + b.m_y;
	real z = a.m_z + b.m_z;
	return { x, y, z };
}

// Returns a Vector resulting from substracting vector B from A.
// A - B
Vector Vector::subVectors(const Vector& a, const Vector& b) {
	real x = a.m_x - b.m_x;
	real y = a.m_y - b.m_y;
	real z = a.m_z - b.m_z;
	return { x, y, z };
}

// Returns the vector result from negating all vector A xyz values.
Vector Vector::negate(const Vector& a) {
	return { -a.m_x, -a.m_y, -a.m_z };
}

// Returns vector resulting from mutiplying vector A by X.
Vector Vector::scalarMultiplication(const Vector& a, const real& x) {
	return { x * a.m_x, x * a.m_y, x * a.m_z };
}

// Returns vector resulting from dividing vector A by X.
Vector Vector::scalarDivision(const Vector& a, const real& x) {
	return { a.m_x / x, a.m_y / x, a.m_z / x };
}

// Returns the magnitude of vector A
// |A|
real Vector::magnitude(const Vector& a) {
	return sqrt(a.m_x * a.m_x + a.m_y * a.m_y + a.m_z * a.m_z);
}

// Returns the vector A normalized.
Vector Vector::normalize(const Vector& a) {
	real x = magnitude(a);
	return { a.m_x / x, a.m_y / x, a.m_z / x };
}

// Returns the dot product of vectors A and B.
// A * B
real Vector::dotProduct(const Vector& a, const Vector& b) {
	return { a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z };
}

// Returns the cross product of vectors A and B.
// A x B
Vector Vector::crossProduct(const Vector& a, const Vector& b) {
	return { a.m_y * b.m_z - a.m_z * b.m_y,
		a.m_z * b.m_x - a.m_x * b.m_z,
		a.m_x * b.m_y - a.m_y * b.m_x };
}

// Returns the reflected vector resulting from the IN vector 
// intersecting an object at a point with NORMAL vector.
Vector Vector::reflected(Vector& in, Vector& normal) {
	return subVectors(in, scalarMultiplication(normal, 2 * dotProduct(in, normal)));
}