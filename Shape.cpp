#include "Shape.h"
#include <iostream>

using namespace rt;

// Default constructor. Assigns a random ID, and an identity tranform matrix.
Shape::Shape() { b_id = rand() % 1000; b_transform = Matrix::identity(4); b_inv_tx = Matrix::inverse(b_transform); }

// Object constructed with specified material and transform matrix. ID is random.
Shape::Shape(Material m, Matrix t) : b_material{ m }, b_transform{ t } { b_id = rand() % 1000; b_inv_tx = Matrix::inverse(b_transform); }

// Destructor.
Shape::~Shape() {}

// Returns an Object ID.
int Shape::getID() const { return b_id; }

void Shape::intersection(const Ray& r, std::vector<real>& hits, std::unordered_map<real, Shape*>& hmap) {
	std::cout << "shape intersecion function\n";
}

void Shape::normal_at(Point& hitpoint, Matrix& inv_tx, Vector& normal) {
	std::cout << "shape normal_at function\n";
	//return (Vector());
}
