/*
* Interface file for the Ray Tracer Matrices.
*
*
* This file contains the definitions for the Matrix class,
* which are used to transform all the objects within a 3D scene.
*
*/

#ifndef MATRIX_H
#define MATRIX_H

#include "Point.h"
#include "Vector.h"
#include "Precision.h"

namespace rt {

	/*
	* The Ray Tracer implements Matrix Transformations to manipulate
	* points and vectors, and scale, translate, rotate, and sheare
	* objects within a 3D scene.
	*/

	class Matrix {
	public:

		// Matrix dimension parameters.
		int m_rows;
		int m_columns;
		int size;

		// A matrix is implemented as a two dimensional array of real numbers.
		real** matrix;
		void create_matrix(int rows, int columns);

		// Default constructor is an identity matrix of size 4.
		Matrix();

		// Construct an square matrix with specified size.
		Matrix(int size);

		// Construct matrix with specified dimensions.
		Matrix(int rows, int columns);

		// Returns the value stored at position 'xy' in the matrix.
		real valueAt(int x, int y);

		// Evaluates if two matrices are identical.
		static bool identical_matrix(const Matrix& a, const Matrix& b);

		// Creates an identity matrix of the specified size.
		static Matrix identity(int size);

		// Multiplies two matrices - AxB - and returns the resulting matrix.
		static Matrix multiplyMatrices(const Matrix& a, const Matrix& b);

		// Multiplies a matrix by a Point.
		static Point multiply_point(const Matrix& a,const Point& b);

		// Multiplies a matrix by a Vector.
		static Vector multiply_vector(const Matrix& a, const Vector& b);

		// Multiplies a matrix by an identity matrix, and returns the resulting matrix.
		static Matrix multiplyIdentity(const Matrix& a);

		// Returns the transpose of matrix A.
		static Matrix transpose(Matrix a);

		// Returns the submatrix XY of matrix M.
		static Matrix submatrix(const Matrix& m, int x, int y);

		// Returns the determinant of matrix M.
		static real determinant(const Matrix& m);

		// Returns the Cofactor Matrix of matrix M.
		static Matrix cofactorMatrix(const Matrix& m);

		// Returns the Inverse matrix of matrix A.
		static Matrix inverse(const Matrix& a);

		// Returns the Translation matrix to point P.
		static Matrix translation_matrix_p(const Point& p);

		// Returns the Translation matrix to vector V.
		static Matrix translation_matrix_v(const Vector& v);

		// Returns an identity matrix scaled by point P xyz values.
		static Matrix scaling_matrix_p(const Point& p);

		// Returns an identity matrix scaled by vector V xyz values.
		static Matrix scaling_matrix_v(const Vector& v);

		// Returns an identity matrix rotated around the X axis by D degrees.
		static Matrix rotation_matrix_x(const real& d);

		// Returns an identity matrix rotated around the Y axis by D degrees.
		static Matrix rotation_matrix_y(const real& d);

		// Returns an identity matrix rotated around the Z axis by D degrees.
		static Matrix rotation_matrix_z(const real& d);

		// Returns an identity matrix sheared according to the specified values.
		static Matrix shearing_matrix(const real& xy, const real& xz, const real& yx,
			const real& yz, const real& zx, const real& zy);

		// Returns a transformation matrix that orients the 3D scene relative to the viewer.
		static Matrix view_transform(Point& from, Point& to,const Vector& up);

	};
}


#endif // !MATRIX_H
