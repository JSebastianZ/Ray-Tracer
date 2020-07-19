#include "Tests.h"

using namespace rt;
/*
void Tests::create_tuple()
{
	Tuple tuple{};
	tuple.m_x = 4.3;
	tuple.m_y = -4.2;
	tuple.m_z = 3.1;
	tuple.m_w = 0.0;
	std::cout << tuple.m_x << ' ' << tuple.m_y << ' ' << tuple.m_z << '\n';
};

void Tests::test_point() {
	Point point{};
	point.m_x = 1;
	point.m_y = 2;
	point.m_z = 3;
	std::cout << point.m_x << ' ' << point.m_y << ' ' << point.m_z << ' ' << point.m_w << '\n';
};

void Tests::test_vector() {
	Vector v{};
	v.m_x = 1;
	v.m_y = 2;
	v.m_z = 3;
	std::cout << v.m_x << ' ' << v.m_y << ' ' << v.m_z << ' ' << v.m_w << '\n';
};

void Tests::vector_norm() {
	Vector v(1, 1, 1);
	Vector normal(Vector::normalize(v));
	std::cout << "vector normal:\n";
	std::cout << normal.m_x << ' ' << normal.m_y << ' ' << normal.m_z << '\n';
}

void Tests::test_color() {
	Color c1{};
	c1.m_r = 1;
	c1.m_g = 0.2;
	c1.m_b = 0.4;
	Color c2{};
	c2.m_r = 0.9;
	c2.m_g = 1;
	c2.m_b = 0.1;

	Color c;
	c.multColors(c1, c2);
	std::cout << c.m_r << ' ' << c.m_g << ' ' << c.m_b << '\n';

};

void Tests::test_canvas() {
	Canvas canvas(10,10);
	Color c(1,1,1);
	std::cout << c.m_r << ' ' << c.m_g << ' ' << c.m_b << '\n';
	canvas.write_pixel(0, 0, c);
	Color b;

	std::cout << b.m_r << ' ' << b.m_g << ' ' << b.m_b << '\n';
	b = canvas.get_pixel(0, 0);
	std::cout << b.m_r << ' ' << b.m_g << ' ' << b.m_b << '\n';
}

void Tests::test_translation() {
	Point p(-3, 4, 5);
	Point b(5, -3, 2);
	Matrix transform(Matrix::translation_matrix_p(b));
	//Point a(Matrix::multiply_point(transform, p));
	Point a = Transform::translate_point(p, b);
	std::cout << a.m_x << ' ' << a.m_y << ' ' << a.m_z << '\n';
}

void Tests::mult_by_inv_transmatrix() {
	Point b(5, -3, 2);
	Point p(-3, 4, 5);
	Matrix m(Matrix::translation_matrix_p(Point(5, -3, 2)));
	Matrix inv(Matrix::inverse(m));
	Point a(Matrix::multiply_point(inv, p));
	std::cout << a.m_x << ' ' << a.m_y << ' ' << a.m_z << '\n';

}

void Tests::scale_point() {
	Point s(2, 3, 4);
	Matrix scaling(Matrix::scaling_matrix_p(s));
	Point p(-4, 6, 8);
	Point transform(Matrix::multiply_point(scaling, p));
	std::cout << transform.m_x << ' ' << transform.m_y << ' ' << transform.m_z << '\n';
}

void Tests::scale_vector() {
	Point s(2, 3, 4);
	Matrix scaling(Matrix::scaling_matrix_p(s));
	Vector v(-4, 6, 8);
	Vector transform(Matrix::multiply_vector(scaling, v));
	std::cout << transform.m_x << ' ' << transform.m_y << ' ' << transform.m_z << '\n';
}

void Tests::mult_by_inv_scalematrix() {
	Point s(2, 3, 4);
	Matrix scaling(Matrix::scaling_matrix_p(s));
	Matrix inv(Matrix::inverse(scaling));
	Vector v(-4, 6, 8);
	Vector transform(Matrix::multiply_vector(inv, v));
	std::cout << transform.m_x << ' ' << transform.m_y << ' ' << transform.m_z << '\n';
}

void Tests::reflection() {
	Matrix transform(Matrix::scaling_matrix_p(Point(-1, 1, 1)));
	Point p(2, 3, 4);
	Point a(Matrix::multiply_point(transform, p));
	std::cout << a.m_x << ' ' << a.m_y << ' ' << a.m_z << '\n';
}

void Tests::rotate_point() {
	Point p(0, 1, 0);
	Matrix half_quarter(Matrix::rotation_matrix_x(pi / 4));
	Matrix full_quarter(Matrix::rotation_matrix_x(pi / 2));
	Point a(Matrix::multiply_point(half_quarter, p));
	Point b(Matrix::multiply_point(full_quarter, p));
	std::cout << "rotate around X axis tests:\n";
	std::cout << a.m_x << ' ' << a.m_y << ' ' << a.m_z << '\n';
	std::cout << b.m_x << ' ' << b.m_y << ' ' << b.m_z << '\n';
	Point w(0, 0, 1);
	Matrix hq(Matrix::rotation_matrix_y(pi / 4));
	Matrix fq(Matrix::rotation_matrix_y(pi / 2));
	Point c = Matrix::multiply_point(hq, w);
	Point d = Matrix::multiply_point(fq, w);
	std::cout << "rotate around Y axis tests:\n";
	std::cout << c.m_x << ' ' << c.m_y << ' ' << c.m_z << '\n';
	std::cout << d.m_x << ' ' << d.m_y << ' ' << d.m_z << '\n';
	Matrix halfq(Matrix::rotation_matrix_z(pi / 4));
	Matrix fullq(Matrix::rotation_matrix_z(pi / 2));
	Point e = Matrix::multiply_point(halfq, p);
	Point f = Matrix::multiply_point(fullq, p);
	std::cout << "rotate around Z axis tests:\n";
	std::cout << e.m_x << ' ' << e.m_y << ' ' << e.m_z << '\n';
	std::cout << f.m_x << ' ' << f.m_y << ' ' << f.m_z << '\n';
}

void Tests::inv_xrot() {
	Point p(0, 1, 0);
	Matrix half_quarter(Matrix::rotation_matrix_x(pi / 4));
	Matrix inv(Matrix::inverse(half_quarter));
	Point transform(Matrix::multiply_point(inv, p));
	std::cout << transform.m_x << ' ' << transform.m_y << ' ' << transform.m_z << '\n';
}

void Tests::shearing_tests() {
	Matrix transform(Matrix::shearing_matrix(0, 1, 0, 0, 0, 0));
	Point p(2, 3, 4);
	Point r(Matrix::multiply_point(transform, p));
	std::cout << r.m_x << ' ' << r.m_y << ' ' << r.m_z << '\n';

}

void Tests::chain_transform() {
	Point p(1, 0, 1);
	Matrix a(Matrix::rotation_matrix_x(pi / 2));
	Matrix b(Matrix::scaling_matrix_p(Point(5, 5, 5)));
	Matrix c(Matrix::translation_matrix_p(Point(10, 5, 7)));
	Point p4(Matrix::multiply_point(c, (Matrix::multiply_point(b, (Matrix::multiply_point(a, p))))));
	std::cout << p4.m_x << ' ' << p4.m_y << ' ' << p4.m_z << '\n';
}

void Tests::material_test() {
	Material m{Material()};
	std::cout << m.m_color.m_r << ' ' << m.m_ambient << '\n';
}

//void Tests::shape_test() {
//	Shape s;
//	std::cout << s.getID() << '\n';
//	std::cout << s.m_id << ' ' << s.m_material.m_ambient << ' ' << s.m_transform.matrix[3][3] << '\n';
//}


void Tests::sphere_test() {
	Sphere s{Sphere()};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << s.b_transform.matrix[i][j];
		}
		std::cout << '\n';
	}

	std::cout << "sphere ID: " << s.b_id << '\n';
	std::cout << "sphere origin: " << s.m_center.m_x << ' ' << s.m_center.m_y << ' ' << s.m_center.m_z << '\n';
	std::cout << "sphere radius: " << s.m_radius << '\n';

}

void Tests::trans_sphere() {
	Sphere s;
	Matrix t(Matrix::translation_matrix_p(Point(2, 3, 4)));
	s.b_transform = t;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << s.b_transform.matrix[i][j];
		}
		std::cout << '\n';
	}
}

void Tests::ray_test() {
	Ray r;
	std::cout << r.m_origin.m_x << ' ' << r.m_origin.m_y << ' ' << r.m_origin.m_z << '\n';
}

void Tests::base_ix() {
	Ray r;
	Sphere s;
	//s.intersection(r);
}

void Tests::sphere_ix() {
	Ray r(Point(0, 0, 5), Vector(0, 0, 1));
	Sphere s;
	s.b_transform = Matrix::translation_matrix_p(Point(3, 4, 5));
	//s.intersection(r);
	//std::cout << g_ix_points[0] << ' ' << g_ix_points[1] << "\n";
}

void Tests::translate_ray() {
	Ray r(Point(1, 2, 3), Vector(0, 1, 0));
	Matrix m(Matrix::translation_matrix_p(Point(3, 4, 5)));
	Matrix s(Matrix::scaling_matrix_p(Point(2, 3, 4)));
	Ray r2(Ray::transform(r, s));
	std::cout << "r2 origin: " << r2.m_origin.m_x << ' ' << r2.m_origin.m_y << ' ' << r2.m_origin.m_z << '\n';
	std::cout << "r2 direction: " << r2.m_direction.m_x << ' ' << r2.m_direction.m_y << ' ' << r2.m_direction.m_z << '\n';
}

void Tests::simple_tracer() {
	Tracer t;
	Ray r(Point(0, 0, -5), Vector(0, 0, 1));
	Sphere s;
	s.intersection(r,t.ix_points,t.ix_shape_map);
	for (auto& element : t.ix_shape_map) {
		std::cout << element.first << " " << element.second <<  "\n";
	}

}

void Tests::light_test() {
	Light l(Point(0, 0, 0), Color(1, 1, 1));
	std::cout << l.m_intensity.m_r << ' ' << l.m_intensity.m_g << ' ' << l.m_intensity.m_b << '\n';
	std::cout << l.m_source.m_x << ' ' << l.m_source.m_y << ' ' << l.m_source.m_z << '\n';
}

void Tests::world_test() {
	Ray r(Point(0, 0, -5), Vector(0, 0, 1));
	Tracer t;
	World w = World();

	for (int i = 0; i < w.m_objects.size(); i++) {
		w.m_objects[i]->intersection(r, t.ix_points, t.ix_shape_map);
	}

	std::cout << "numero de intersecciones: " << t.ix_points.size() << '\n';

	for (int i = 0; i < t.ix_points.size(); i++) {
		std::cout << t.ix_points[i] << '\n';
	}

	std::cout << "mapeo\n";
	for (auto& t : t.ix_shape_map) {
		std::cout << t.first << " " << t.second->b_id << "\n";

	}
};

void Tests::color_at() {
	Ray r(Point(0, 0, -5), Vector(0, 0, 1));
	Tracer t;
	World w = World();
	Color color;
	color = t.color_at(w, r);
	std::cout << "color at \n";
	std::cout << color.m_r << ' ' << color.m_g << ' ' << color.m_b << '\n';
}

void Tests::view_tx() {
	Point from(1, 3, 2);
	Point to(4, -2, 8);
	Vector up(1, 1, 0);
	Matrix t(Matrix::view_transform(from, to, up));
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			std::cout << t.matrix[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

void Tests::icamera() {
	Camera c(160, 120,pi/2);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << c.m_transform.matrix[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

void Tests::pixel_size() {
	Camera c(125, 200, pi / 2);
	std::cout << "pixel size: " << c.pixel_size() << '\n';
}

void Tests::camera_r4p() {
	Camera c(201, 101, pi / 2);
	c.m_transform = Matrix::multiplyMatrices(Matrix::rotation_matrix_y(pi / 4),
		Matrix::translation_matrix_p(Point(0, -2, 5)));
	Ray r = c.ray_for_pixel(100, 50);
	std::cout << "ray origin: " << r.m_origin.m_x << ',' << r.m_origin.m_y << ',' << r.m_origin.m_z << '\n';
	std::cout << "ray direction: " << r.m_direction.m_x << ',' << r.m_direction.m_y << ',' << r.m_direction.m_z << '\n';
}

void Tests::render() {
	World w = World();
	Camera c(11, 11, pi / 2);
	Point from(0, 0, -5);
	Point to(0, 0, 0);
	Vector up(0, 1, 0);
	c.m_transform = Matrix::view_transform(from, to, up);
	Tracer t;
	t.render(c, w);
	
}

void Tests::pit7() {
	World w = World();
	Tracer t;
	Camera camera = Camera(128, 128, (pi/3));
	Point from = Point(0, 1.5, -5);
	Point to = Point(0, 1, 0);
	Vector up = Vector(0, 1, 0);
	camera.m_transform = Matrix::view_transform(from, to, up);
	camera.m_inv_tx = Matrix::inverse(camera.m_transform);
	t.render(camera, w);
}

void Tests::verga() {
	Color black{ 0,0,0 };
	Point p{ 0,0,140 };
	Color ctest(1,0.2,1);
	Material material(ctest, 0.1, 0.9, 0.9, 200.0);
	Color light_color(Color(1, 1, 1));
	Point light_position(Point(-50, 50, -10));
	Light light(light_position, light_color);
	Canvas knvas = Canvas(128, 128);
	int knvas_depth = 128;
	Point scale{ 20 , 20, 20};
	Matrix t(Matrix::scaling_matrix_p(scale));
	Matrix transform(Matrix::multiplyMatrices(Matrix::translation_matrix_p(p), t));
	Matrix trans(Matrix::translation_matrix_p(p));
	Sphere sphere = Sphere(2*scale.m_x, p, material, transform);
	Tracer tracer;
	for (int x = (-knvas.m_width / 2), i = 0; x < (knvas.m_width / 2); x++, i++) {
		for (int y = (-knvas.m_height / 2), j = 0; y < (knvas.m_height / 2); y++, j++) {
			// Creo el rayo y determino si le pega a la esfera o no
			Ray ray(Point(0, 0, 0), (Vector::normalize(Point::subPoints(Point(x, y, knvas_depth), Point(0, 0, 0)))));
			Ray ray2(ray.m_origin, ray.m_direction);
			Ray r = Ray::transform(ray2, Matrix::inverse(sphere.b_transform));
			sphere.intersection(r, tracer.ix_points, tracer.ix_shape_map);
			tracer.sort_ix(tracer.ix_points);
			tracer.m_hit3D = tracer.min_hit_point(tracer.ix_points);
			tracer.ix_points.clear();
			// Defino el vector EYE como el negativo del rayo
			tracer.eye = (Vector::negate(r.m_direction));
			// Ahora necesito los PUNTOS de INTERSECCION en 3D Phit = ray.orig + ray.dir * t
			tracer.hit3D = (Ray::position(r, tracer.m_hit3D));
			Matrix mti(Matrix::inverse(transform));
			// Calculo del vector Normal al hit point choke
			Point object_point(Matrix::multiply_point(mti, tracer.hit3D));
			Vector object_normal(Point::subPoints(object_point, Point(0, 0, 0)));
			Vector world_normal(Matrix::multiply_vector(Matrix::transpose(Matrix::inverse(transform)), object_normal));
			world_normal.m_w = 0;
			tracer.normal = (Vector::normalize(world_normal));
			// Lightning
			Color c(tracer.lighting(sphere.b_material, light, false));
			c.m_r = static_cast<int>(255.999 * c.m_r);
			if (c.m_r > 255) c.m_r = 255;
			c.m_g = static_cast<int>(255.999 * c.m_g);
			if (c.m_g > 255) c.m_g = 255;
			c.m_b = static_cast<int>(255.999 * c.m_b);
			if (c.m_b > 255) c.m_b = 255;
			if (tracer.m_hit3D >= 0) {
				knvas.write_pixel(i, j, c);
			}
			else {
				knvas.write_pixel(i, j, black);
			}
		}
		tracer.m_hit3D = 0;
	}
	std::ofstream tcfile;
	tcfile.open("maricon.ppm");
	tcfile << "P3\n" << knvas.m_width << ' ' << knvas.m_height << '\n' << 255 << '\n';
	for (int j = knvas.m_height - 1; j >= 0; --j) {
		for (int i = 0; i < knvas.m_width; ++i) {
			tcfile << knvas.get_pixel(i, j).m_r << ' ' << knvas.get_pixel(i, j).m_g << ' ' << knvas.get_pixel(i, j).m_b << '\n';
		}
	}
	tcfile.close();
	std::cout << "maricon" << '\n';
}

void Tests::circle() {
	real canvas_pixels = 100;
	Canvas canvas(canvas_pixels, canvas_pixels);
	Color color(1, 0, 0);
	Color black(0, 0, 0);
	real wallz = 10;
	real wallsize = 7;
	real half = wallsize / 2;
	real pixel_size = wallsize / canvas_pixels;
	Sphere sphere;
	Tracer t;
	for (int y = 0; y < canvas_pixels; y++) {
		real worldy = half - pixel_size * y;
		for (int x = 0; x < canvas_pixels; x++) {
			real worldx = -half + pixel_size * x;
			Point position(worldx,worldy,wallz);
			Ray r(Point(0, 0, -5), Vector::normalize(Point::subPoints(position, Point(0, 0, -5))));
			sphere.intersection(r, t.ix_points, t.ix_shape_map);
			std::cout << t.ix_points.size() << '\n';
			if (t.ix_points[0] != infinite) {
				color.m_r = static_cast<int>(255.999 * color.m_r);
				if (color.m_r > 255) color.m_r = 255;
				color.m_g = static_cast<int>(255.999 * color.m_g);
				if (color.m_g > 255) color.m_g = 255;
				color.m_b = static_cast<int>(255.999 * color.m_b);
				if (color.m_b > 255) color.m_b = 255;
				canvas.write_pixel(x, y, color);
				t.ix_points.clear();
			}
			else {
				canvas.write_pixel(x, y, black);
				t.ix_points.clear();
			}
		}
	}

	std::ofstream tcfile;
	tcfile.open("maricon.ppm");
	tcfile << "P3\n" << canvas.m_width << ' ' << canvas.m_height << '\n' << 255 << '\n';
	for (int j = canvas.m_height - 1; j >= 0; --j) {
		for (int i = 0; i < canvas.m_width; ++i) {
			tcfile << canvas.get_pixel(i, j).m_r << ' ' << canvas.get_pixel(i, j).m_g << ' ' << canvas.get_pixel(i, j).m_b << '\n';
		}
	}
	tcfile.close();
	std::cout << "maricon" << '\n';
}

void Tests::scale_sphere_ix() {
	Tracer t;
	Ray r(Point(0, 0, -5), Vector(0, 0, 1));
	Sphere s;
	s.b_transform = Matrix::translation_matrix_p(Point(5,0,0));
	s.intersection(r, t.ix_points, t.ix_shape_map);
	std::cout << t.ix_points.size() << '\n';
}

void Tests::pit7_sincamera() {
	World w;
	Tracer t;
	std::ofstream tcfile;
	tcfile.open("test.ppm");
	tcfile << "P3\n" << 150 << ' ' << 150 << '\n' << 255 << '\n';
	Canvas image = Canvas(150, 150);
	for (int x = 0; x < 150; x++) {
		for (int y = 0; y < 150; y++) {
			Ray ray (Point(0, 0, 0), (Vector::normalize(Point::subPoints(Point(x, y, 150), Point(0, 0, 0)))));
			Color color = t.color_at(w, ray);
			t.ix_points.clear();				// if not clear the tracer array keeps adding records.
			t.ix_shape_map.clear();
			image.write_pixel(x, y, color);
			real r = color.m_r * 255;
			if (r > 255) r = 255;
			real g = color.m_g * 255;
			if (g > 255) g = 255;
			real b = color.m_b * 255;
			if (b > 255) b = 255;
			tcfile << r << ' ' << g << ' ' << b << '\n';
		}
	}
	tcfile.close();
	std::cout << "maricon" << '\n';
}

void Tests::plane_const_normal() {
	Tracer t;
	Plane plane;
	Point p(-5, 0, 150);
	plane.normal_at(p, plane.b_transform, t.normal);
	std::cout << t.normal.m_x << ',' << t.normal.m_y << ',' << t.normal.m_z << '\n';
}

void Tests::ray_plane_ix_parallel_or_coplanar() {
	Tracer t;
	Plane plane;
	Ray r(Point(0, 10, 0), Vector(0, 0, 1));
	plane.intersection(r, t.ix_points, t.ix_shape_map);
	std::cout << "el tamanho del array ixs es: " << t.ix_points.size() << '\n';
	std::cout << "los hits estan en: " << '\n';
	for (int i = 0; i < t.ix_points.size(); i++) {
		std::cout << t.ix_points[i] << '\n';
	}
}

void Tests::ix_plane_above() {
	Tracer t;
	Plane* p = new Plane();
	Ray r(Point(0, 1, 0), Vector(0, -1, 0));
	p->intersection(r, t.ix_points, t.ix_shape_map);
	std::cout << "el tamanho del array ixs es: " << t.ix_points.size() << '\n';
	std::cout << "los hits estan en: " << '\n';
	for (int i = 0; i < t.ix_points.size(); i++) {
		std::cout << t.ix_points[i] << '\n';
	}
}

void Tests::pit9() {

	World w = World();
	Tracer t;
	Camera camera = Camera(256, 256, (pi / 3));
	Point from = Point(0, 1.5, -5);
	Point to = Point(0, 1, 0);
	Vector up = Vector(0, 1, 0);
	camera.m_transform = Matrix::view_transform(from, to, up);
	camera.m_inv_tx = Matrix::inverse(camera.m_transform);
	t.render(camera, w);
}

void Tests::stripe_pattern() {
	Color white(1, 1, 1);
	Color black(0, 0, 0);
	Color test;
	Material m;
	//m.m_stripe_pattern.push_back(white);
	//m.m_stripe_pattern.push_back(black);

	Point p(-1.1, 0, 0);
	//Color verga = Color::stripe_at(m.m_stripe_pattern, p);
	//std::cout << verga.m_r << ',' << verga.m_g << ',' << verga.m_b << '\n';
}

void Tests::lighting_pattern() {
	Tracer t;
	t.hit3D = Point(0.9, 0, 0);
	Material m;
	Color white(1, 1, 1);
	Color black(0, 0, 0);
	//m.m_stripe_pattern.push_back(white);
	//m.m_stripe_pattern.push_back(black);
	//m.m_pattern = true;
	m.m_ambient = 1;
	m.m_diffuse = 0;
	m.m_specular = 0;
	Vector eyev(0, 0, -1);
	Vector normalv(0, 0, -1);
	Light light(Point(0, 0, -10), Color(1, 1, 1));
	Color c1, c2;
	c1 = t.lighting(m, light, false);
	std::cout << "c1: " << c1.m_r << ',' << c1.m_g << ',' << c1.m_b << '\n';
	t.hit3D = Point(1.1, 0, 0);
	c2 = t.lighting(m, light, false);
	std::cout << "c2: " << c2.m_r << ',' << c2.m_g << ',' << c2.m_b << '\n';
}
*/

/*
void Tests::reflection() {

	Plane shape;
	Ray r(Point(0, 1, -1), Vector(0, -0.7071678, 0.7071678));
	Tracer t;
	shape.intersection(r, t.ix_points, t.ix_shape_map);
	for (int i = 0; i < t.ix_points.size(); i++) {
		std::cout << t.ix_points[i] << '\n';
	}
	Light l(Point(-10, 10, -10), Color(1, 1, 1));
	Vector reflex = Vector::reflected(r.m_direction, shape.m_plane_normal);
	std::cout << reflex.m_x << ',' << reflex.m_y << ',' << reflex.m_z << '\n';
}
*/

void Tests::reflux() {
	World w = World();
	Tracer t;
	Camera camera = Camera(512, 512, (pi/1.5));
	Point from = Point(-1, 1.5, -5);
	Point to = Point(0, 1, 0);
	Vector up = Vector(0, 1, 0);
	camera.m_transform = Matrix::view_transform(from, to, up);
	camera.m_inv_tx = Matrix::inverse(camera.m_transform);
	t.render(camera, w);
}