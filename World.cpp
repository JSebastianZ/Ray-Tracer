#include "World.h"
#include <iostream>

using namespace rt;

World::World() {
	Color white(1, 1, 1);
	Color black(0, 0, 0);

	// Default world
	/*
	m_light_source = Light(Point(-10, 10, -10), Color(1, 1, 1));
	Sphere* s1 = new Sphere();
	s1->m_radius = 1;
	s1->b_material = Material();
	s1->b_material.m_color = Color(1,0.2,1);
	s1->b_material.m_diffuse = 0.7;
	s1->b_material.m_specular = 0.2;
	s1->b_material.m_shininess = 200;
	s1->b_transform = Matrix::translation_matrix_p(Point(0, 0, 0));


	Sphere* s2 = new Sphere();
	s2->m_diameter = 1;
	s2->m_radius = 0.5;
	s2->b_material = Material();
	s2->b_material.m_ambient = 0.1;
	s2->b_material.m_diffuse = 0.9;
	s2->b_material.m_specular = 0.9;
	s2->b_material.m_shininess = 200;
	s2->b_material.m_color = Color(1, 1, 1);
	s2->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(0, 0, 0)),
		Matrix::scaling_matrix_p(Point(0.5, 0.5, 0.5)));
	m_objects.push_back(s1);
	//std::cout << s1->b_material.m_color.m_r << ',' << s1->b_material.m_color.m_g << ',' << s1->b_material.m_color.m_b << '\n';
	//m_objects.push_back(s2);
	//std::cout << s2->b_material.m_color.m_r << ',' << s2->b_material.m_color.m_g << ',' << s2->b_material.m_color.m_b << '\n';
	*/

	// World for reflection
	m_light_source = Light(Point(-10, 10, -10), Color(0.5, 0.5, 0.5));
	Plane* p1 = new Plane();
	Plane* p2 = new Plane();
	Plane* p3 = new Plane();
	Sphere* s4 = new Sphere();
	Sphere* s5 = new Sphere();
	Sphere* s6 = new Sphere();

	p1->b_transform = Matrix::translation_matrix_p(Point(0, 0, 0));
	p1->b_inv_tx = Matrix::inverse(p1->b_transform);
	p1->b_material.m_color = Color(1, 1, 1);
	p1->b_material.m_specular = 0;
	p1->b_material.m_reflective = 1;
	p1->b_material.m_pattern.m_pattern_exist = true;
	p1->b_material.m_pattern.pt = CHECKER;
	p1->b_material.m_pattern.m_stripe_pattern.push_back(white);
	p1->b_material.m_pattern.m_stripe_pattern.push_back(black);
	p1->b_material.m_pattern.b_pattern_tx = Matrix::scaling_matrix_p(Point(2, 2, 2));
	p1->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(p1->b_material.m_pattern.b_pattern_tx);

	p2->b_transform = //Matrix::translation_matrix_p(Point(0, 0, 0));
		Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(2, 0, 8)),
			(Matrix::multiplyMatrices(Matrix::rotation_matrix_y(pi/4),
				Matrix::rotation_matrix_x(pi / 2)))
		);
	p2->b_inv_tx = Matrix::inverse(p2->b_transform);
	p2->b_material.m_color = Color(1, 1, 1);
	p2->b_material.m_specular = 0;
	p2->b_material.m_reflective = 1;
	p2->b_material.m_ambient = 0;
	p2->b_material.m_shininess = 0;
	p2->b_material.m_diffuse = 1;
	p2->b_material.m_pattern.m_pattern_exist = false;
	p2->b_material.m_pattern.pt = STRIPE;
	p2->b_material.m_pattern.m_stripe_pattern.push_back(white);
	p2->b_material.m_pattern.m_stripe_pattern.push_back(black);
	p2->b_material.m_pattern.b_pattern_tx = Matrix::scaling_matrix_p(Point(2, 2, 2));
	p2->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(p2->b_material.m_pattern.b_pattern_tx);

	p3->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(2, 0, 8)),
		(Matrix::multiplyMatrices(Matrix::rotation_matrix_y(-pi / 4),
			Matrix::rotation_matrix_x(pi / 2)))
	);
	p3->b_inv_tx = Matrix::inverse(p3->b_transform);
	p3->b_material.m_color = Color(1, 1, 1);
	p3->b_material.m_specular = 0;
	p3->b_material.m_reflective = 1;
	p3->b_material.m_pattern.m_pattern_exist = false;
	p3->b_material.m_pattern.pt = STRIPE;
	p3->b_material.m_pattern.m_stripe_pattern.push_back(white);
	p3->b_material.m_pattern.m_stripe_pattern.push_back(black);
	p3->b_material.m_pattern.b_pattern_tx =	Matrix::scaling_matrix_p(Point(2, 2, 2));
	p3->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(p3->b_material.m_pattern.b_pattern_tx);

	s4->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(-0.15, 1, 0.5)),
		Matrix::rotation_matrix_z(pi / 2));
	s4->b_inv_tx = Matrix::inverse(s4->b_transform);
	s4->b_material = Material();
	s4->b_material.m_color = Color(0.1, 1, 0.5);
	s4->b_material.m_diffuse = 0.7;
	s4->b_material.m_specular = 0.3;
	s4->b_material.m_reflective = 0;
	s4->b_material.m_pattern.m_pattern_exist = false;
	//s4->b_material.m_pattern.pt = RING;
	//s4->b_material.m_pattern.m_stripe_pattern.push_back(white);
	//s4->b_material.m_pattern.m_stripe_pattern.push_back(black);
	//s4->b_material.m_pattern.b_pattern_tx = Matrix::multiplyMatrices(
	//	Matrix::scaling_matrix_p(Point(0.2, 0.2, 0.2)),
	//	Matrix::rotation_matrix_z(pi / 2)
	//);
	//s4->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s4->b_material.m_pattern.b_pattern_tx);

	s5->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(1.5, 0.5, -0.5)),
		Matrix::scaling_matrix_p(Point(0.5, 0.5, 0.5)));
	s5->b_inv_tx = Matrix::inverse(s5->b_transform);
	s5->b_material = Material();
	s5->b_material.m_color = Color(0.5, 1, 0.1);
	s5->b_material.m_diffuse = 0.7;
	s5->b_material.m_specular = 0.3;
	s5->b_material.m_reflective = 0;
	s5->b_material.m_pattern.m_pattern_exist = false;
	s5->b_material.m_pattern.pt = CHECKER;
	s5->b_material.m_pattern.m_stripe_pattern.push_back(white);
	s5->b_material.m_pattern.m_stripe_pattern.push_back(black);
	s5->b_material.m_pattern.b_pattern_tx = Matrix::translation_matrix_p(Point(1, 1, 1));
	s5->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s5->b_material.m_pattern.b_pattern_tx);

	s6->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(-1.5, 0.33, -0.75)),
		Matrix::scaling_matrix_p(Point(0.33, 0.33, 0.33)));
	s6->b_inv_tx = Matrix::inverse(s6->b_transform);
	s6->b_material = Material();
	s6->b_material.m_color = Color(1, 0.8, 0.1);
	s6->b_material.m_diffuse = 0.7;
	s6->b_material.m_specular = 0.3;
	s6->b_material.m_reflective = 0;
	s6->b_material.m_pattern.m_pattern_exist = false;
	s6->b_material.m_pattern.pt = GRADIENT;
	s6->b_material.m_pattern.m_stripe_pattern.push_back(white);
	s6->b_material.m_pattern.m_stripe_pattern.push_back(black);
	s6->b_material.m_pattern.b_pattern_tx = Matrix::translation_matrix_p(Point(1, 1, 1));
	s6->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s6->b_material.m_pattern.b_pattern_tx);

	m_objects.push_back(p1);
	m_objects.push_back(p2);
	m_objects.push_back(p3);
	m_objects.push_back(s4);
	m_objects.push_back(s5);
	m_objects.push_back(s6);

	/*	// Default world for PIT_9 test
	m_light_source = Light(Point(0, 1.50, -1), Color(1, 1, 1));
	Plane* p1 = new Plane();
	Plane* p2 = new Plane();
	Sphere* s4 = new Sphere();
	Sphere* s5 = new Sphere();
	Sphere* s6 = new Sphere();

	p1->b_transform = Matrix::translation_matrix_p(Point(0, 0, 0));
	p1->b_inv_tx = Matrix::inverse(p1->b_transform);
	p1->b_material.m_color = Color(1, 0.9, 0.9);
	p1->b_material.m_specular = 0;
	p1->b_material.m_pattern.m_pattern_exist = true;
	p1->b_material.m_pattern.pt = CHECKER;
	p1->b_material.m_pattern.m_stripe_pattern.push_back(white);
	p1->b_material.m_pattern.m_stripe_pattern.push_back(black);
	p1->b_material.m_pattern.b_pattern_tx = Matrix::scaling_matrix_p(Point(2, 2, 2));
	p1->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(p1->b_material.m_pattern.b_pattern_tx);

	p2->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(0, 0, 5)),
		Matrix::rotation_matrix_x(pi/2));
	p2->b_inv_tx = Matrix::inverse(p2->b_transform);
	p2->b_material.m_color = Color(0, 0, 0);
	p2->b_material.m_specular = 0;
	p2->b_material.m_pattern.m_pattern_exist = true;
	p2->b_material.m_pattern.pt = STRIPE;
	p2->b_material.m_pattern.m_stripe_pattern.push_back(white);
	p2->b_material.m_pattern.m_stripe_pattern.push_back(black);
	p2->b_material.m_pattern.b_pattern_tx = Matrix::scaling_matrix_p(Point(1, 1, 1));
	p2->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(p2->b_material.m_pattern.b_pattern_tx);

	s4->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(-0.15, 1, 0.5)),
		Matrix::rotation_matrix_z(pi/2));
	s4->b_inv_tx = Matrix::inverse(s4->b_transform);
	s4->b_material = Material();
	s4->b_material.m_color = Color(0.1, 1, 0.5);
	s4->b_material.m_diffuse = 0.7;
	s4->b_material.m_specular = 0.3;
	s4->b_material.m_pattern.m_pattern_exist = true;
	s4->b_material.m_pattern.pt = RING;
	s4->b_material.m_pattern.m_stripe_pattern.push_back(white);
	s4->b_material.m_pattern.m_stripe_pattern.push_back(black);
	s4->b_material.m_pattern.b_pattern_tx = Matrix::multiplyMatrices(
		Matrix::scaling_matrix_p(Point(0.2, 0.2, 0.2)),
		Matrix::rotation_matrix_z(pi/2)
		);
	s4->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s4->b_material.m_pattern.b_pattern_tx);

	s5->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(1.5, 0.5, -0.5)),
		Matrix::scaling_matrix_p(Point(0.5, 0.5, 0.5)));
	s5->b_inv_tx = Matrix::inverse(s5->b_transform);
	s5->b_material = Material();
	s5->b_material.m_color = Color(0.5, 1, 0.1);
	s5->b_material.m_diffuse = 0.7;
	s5->b_material.m_specular = 0.3;
	s5->b_material.m_pattern.m_pattern_exist = true;
	s5->b_material.m_pattern.pt = CHECKER;
	s5->b_material.m_pattern.m_stripe_pattern.push_back(white);
	s5->b_material.m_pattern.m_stripe_pattern.push_back(black);
	s5->b_material.m_pattern.b_pattern_tx = Matrix::translation_matrix_p(Point(1, 1, 1));
	s5->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s5->b_material.m_pattern.b_pattern_tx);

	s6->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(-1.5, 0.33, -0.75)),
		Matrix::scaling_matrix_p(Point(0.33, 0.33, 0.33)));
	s6->b_inv_tx = Matrix::inverse(s6->b_transform);
	s6->b_material = Material();
	s6->b_material.m_color = Color(1, 0.8, 0.1);
	s6->b_material.m_diffuse = 0.7;
	s6->b_material.m_specular = 0.3;
	s6->b_material.m_pattern.m_pattern_exist = true;
	s6->b_material.m_pattern.pt = GRADIENT;
	s6->b_material.m_pattern.m_stripe_pattern.push_back(white);
	s6->b_material.m_pattern.m_stripe_pattern.push_back(black);
	s6->b_material.m_pattern.b_pattern_tx = Matrix::translation_matrix_p(Point(1, 1, 1));
	s6->b_material.m_pattern.b_inv_pattern_tx = Matrix::inverse(s6->b_material.m_pattern.b_pattern_tx);

	m_objects.push_back(p1);
	m_objects.push_back(s4);
	m_objects.push_back(s5);
	m_objects.push_back(s6);
	*/
	/*	// Default world for PIT_7

	m_light_source = Light(Point(-10, 10, -10), Color(1, 1, 1));
	Sphere* s1 = new Sphere();
	Sphere* s2 = new Sphere();
	Sphere* s3 = new Sphere();
	Sphere* s4 = new Sphere();
	Sphere* s5 = new Sphere();
	Sphere* s6 = new Sphere();

	s1->m_radius = 1;
	s1->m_diameter = 2;
	s1->b_transform = Matrix::scaling_matrix_p(Point(10, 0.01, 10));
	s1->b_inv_tx = Matrix::inverse(s1->b_transform);
	s1->b_material = Material();
	s1->b_material.m_color = Color(1, 0.9, 0.9);
	s1->b_material.m_specular = 0;

	s2->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(0, 0, 5)),
		Matrix::multiplyMatrices(Matrix::rotation_matrix_y((-1) * pi / 4), 
			Matrix::multiplyMatrices(Matrix::rotation_matrix_x(pi / 2),
				Matrix::scaling_matrix_p(Point(10, 0.01, 10)))));
	s2->b_inv_tx = Matrix::inverse(s2->b_transform);
	s2->b_material = Material();
	s2->b_material.m_color = Color(1, 0.9, 0.9);
	s2->b_material.m_specular = 0;

	s3->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(0, 0, 5)),
		Matrix::multiplyMatrices(Matrix::rotation_matrix_y(pi / 4),
			Matrix::multiplyMatrices(Matrix::rotation_matrix_x(pi / 2),
				Matrix::scaling_matrix_p(Point(10, 0.01, 10)))));
	s3->b_inv_tx = Matrix::inverse(s3->b_transform);
	s3->b_material = Material();
	s3->b_material.m_color = Color(1, 0.9, 0.9);
	s3->b_material.m_specular = 0;

	s4->b_transform = Matrix::translation_matrix_p(Point(-0.5,1,0.5));
	s4->b_inv_tx = Matrix::inverse(s4->b_transform);
	s4->b_material = Material();
	s4->b_material.m_color = Color(0.1, 1, 0.5);
	s4->b_material.m_diffuse = 0.7;
	s4->b_material.m_specular = 0.3;

	s5->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(1.5, 0.5, -0.5)),
		Matrix::scaling_matrix_p(Point(0.5, 0.5, 0.5)));
	s5->b_inv_tx = Matrix::inverse(s5->b_transform);
	s5->b_material = Material();
	s5->b_material.m_color = Color(0.5, 1, 0.1);
	s5->b_material.m_diffuse = 0.7;
	s5->b_material.m_specular = 0.3;

	s6->b_transform = Matrix::multiplyMatrices(Matrix::translation_matrix_p(Point(-1.5, 0.33, -0.75)),
		Matrix::scaling_matrix_p(Point(0.33, 0.33, 0.33)));
	s6->b_inv_tx = Matrix::inverse(s6->b_transform);
	s6->b_material = Material();
	s6->b_material.m_color = Color(1, 0.8, 0.1);
	s6->b_material.m_diffuse = 0.7;
	s6->b_material.m_specular = 0.3;

	m_objects.push_back(s1);
	m_objects.push_back(s2);
	m_objects.push_back(s3);
	m_objects.push_back(s4);
	m_objects.push_back(s5);
	m_objects.push_back(s6);

	*/

}

Sphere World::glass_s() {
	Sphere s = Sphere();
	s.b_material.m_transparecy = 1.0;
	s.b_material.m_refractive_index = 1.5;
	return s;
}