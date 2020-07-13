#include "Transform.h"

using namespace rt;

	/*
	 * --------------------------------------------------------------------------
	 * FUNCTIONS DECLARED IN HEADER:
	 * --------------------------------------------------------------------------
	 */

	Point Transform::translate_point(const Point& p, const Point& b) {
		Matrix m = Matrix::translation_matrix_p(p);
		Point k = Matrix::multiply_point(m, b);
		return k;
	}

	Vector Transform::translate_vector(const Point& p, const Vector& v) {
		Matrix m = Matrix::translation_matrix_p(p);
		Vector k = Matrix::multiply_vector(m, v);
		k.m_w = 0;
		return k;
	}

	Point Transform::scale_point(const Point& p, const Point& b) {
		Matrix m = Matrix::scaling_matrix_p(p);
		Point k = Matrix::multiply_point(m, b);
		return k;
	}

	Vector Transform::scale_vector(const Point& p, const Vector& v) {
		Matrix m = Matrix::scaling_matrix_p(p);
		Vector k = Matrix::multiply_vector(m, v);
		k.m_w = 0;
		return k;
	}

	Point Transform::reflect_point_x(const Point& p) {
		Point n{ -1,1,1 };
		return scale_point(n, p);
	}

	Point Transform::reflect_point_y(const Point& p) {
		Point n{ 1,-1,1 };
		return scale_point(n, p);
	}

	Point Transform::reflect_point_z(const Point& p) {
		Point n{ 1,1,-1 };
		return scale_point(n, p);
	}

	Point Transform::mirror_point(const Point& p) {
		Point n{ -1,-1,-1 };
		return scale_point(n, p);
	}

	Vector Transform::reflect_vector_x(const Vector& v) {
		Point n{ -1,1,1 };
		return scale_vector(n, v);
	}

	Vector Transform::reflect_vector_y(const Vector& v) {
		Point n{ 1,-1,1 };
		return scale_vector(n, v);
	}

	Vector Transform::reflect_vector_z(const Vector& v) {
		Point n{ 1,1,-1 };
		return scale_vector(n, v);
	}

	Vector Transform::mirror_vector(const Vector& v) {
		Point n{ -1,-1,-1 };
		return scale_vector(n, v);
	}

	Point Transform::rotate_point_x(const Point& p, real d) {
		Point r = Matrix::multiply_point(Matrix::rotation_matrix_x(d), p);
		return r;
	}

	Point Transform::rotate_point_y(const Point& p, real d) {
		Point r = Matrix::multiply_point(Matrix::rotation_matrix_y(d), p);
		return r;
	}

	Point Transform::rotate_point_z(const Point& p, real d) {
		Point r = Matrix::multiply_point(Matrix::rotation_matrix_z(d), p);
		return r;
	}

	Vector Transform::rotate_vector_x(const Vector& v, real d) {
		Vector r = Matrix::multiply_vector(Matrix::rotation_matrix_x(d), v);
		r.m_w = 0;
		return r;
	}

	Vector Transform::rotate_vector_y(const Vector& v, real d) {
		Vector r = Matrix::multiply_vector(Matrix::rotation_matrix_y(d), v);
		r.m_w = 0;
		return r;
	}

	Vector Transform::rotate_vector_z(const Vector& v, real d) {
		Vector r = Matrix::multiply_vector(Matrix::rotation_matrix_z(d), v);
		r.m_w = 0;
		return r;
	}

	Point Transform::shearing(const Point& p, const real& xy, const real& xz, const real& yx,
		const real& yz, const real& zx, const real& zy) {
		Point s = Matrix::multiply_point(Matrix::shearing_matrix(xy, xz, yx, yz, zx, zy), p);
		return s;
	}
