#include "Tracer.h"
#include <iostream>
#pragma omp parallel

using namespace rt;

Tracer::Tracer() {}

void Tracer::intersect_world(const World& w,const  Ray& ray) {
	Ray ray2(ray.m_origin, ray.m_direction);					// Ray must be transformed by inverse of shape transform matrix.
	for (int i = 0; i < w.m_objects.size(); i++) {
		Ray r = Ray::transform(ray2, w.m_objects[i]->b_inv_tx);
		//Ray r = Ray::transform(ray2, Matrix::inverse(w.m_objects[i]->b_transform));
		w.m_objects[i]->intersection(r, this->ix_points, this->ix_shape_map);
	}
}

void Tracer::sort_ix(std::vector<real>& ixs) {		
	std::sort(ixs.begin(), ixs.end());
}

real Tracer::min_hit_point(std::vector<real>& ixs) {
	real minimum{ high_number };					 // Initialize "minimum" to a high value for comparison purposes.
	for (int q = 0; q < ixs.size(); q++) {
		if ((ixs[q] < minimum) && (ixs[q] >= 0)) {
			minimum = ixs[q];
		}
		else continue;
	}
	return minimum;
}

void Tracer::prepare_computations( Ray& ray, World& w) {
	sort_ix(this->ix_points);
	this->m_hit3D = min_hit_point(this->ix_points);
	if (this->m_hit3D == high_number) return;
	std::unordered_map<real, Shape*>::const_iterator mapper = ix_shape_map.find(this->m_hit3D);	// Stores object mapped to the hit location
	this->s = mapper->second;															// according to hash map data structure.
	Matrix transform(this->s->b_inv_tx);				// Get the transform matrix of hit object.
	Vector eye = Vector::negate(ray.m_direction);		// Assigns eye vector the negative of ray direction.

	Point hitpoint_3d(Ray::position(ray, this->m_hit3D));		// Computes 3D coordinate of hit point as hitPoint = ray.orig + ray.dir * hit_location.

	s->normal_at(hitpoint_3d,this->s->b_inv_tx,this->normal);	// Normal at hitpoint after shape is transformed.
	///*
	if (Vector::dotProduct(normal, eye) < 0) {			// If normal and eye vector dot product is negative
		ray.m_inside = true;								// then the "viewer" is inside the object.
		normal = Vector::negate(normal);
	}
	else {
		ray.m_inside = false;
	}
	//*/
	Point over_point = Point::addPoints(hitpoint_3d,			// compensate for shadows computation.
		Point::scalarMultiplication
		(Point(normal.m_x, normal.m_y, normal.m_z), epsilon));
	hitpoint_3d = over_point;
	this->hit3D = hitpoint_3d;									// Set hit point, eye, and normal vectors
	this->eye = eye;
	this->normal = normal;
}

bool Tracer::shadowed(World& w, Point& p) {
	Vector v = Point::subPoints(w.m_light_source.m_source, p);
	real distance = Vector::magnitude(v);
	Vector direction = Vector::normalize(v);
	Ray r = Ray(p, direction);
	intersect_world(w, r);
	m_hit3D = min_hit_point(this->ix_points);
	if (m_hit3D < 0 || m_hit3D == high_number) {
		return false;
	}
	else {
		if (m_hit3D < distance) return true; else return false;
	}
}

Color Tracer::lighting(Material& m, Light& light, bool shadow) {
	Color diffuse(0, 0, 0);																// Diffuse reflection component.
	Color specular(0, 0, 0);															// Specular reflection component.
	Color effective_color;
	if (this->s->b_material.m_pattern.m_pattern_exist == true) {						// Compute color based on stripes if defined.
		effective_color = this->s->b_material.m_pattern.stripe_at_object(this->hit3D,this->s->b_inv_tx);
	}
	else {
		effective_color = Color::multColors(m.m_color, light.m_intensity); 			// Combine object surface and light source colors.
	}
	Vector lightv(Vector::normalize(Point::subPoints(light.m_source, this->hit3D))); 			// Vector represents direction of light source.
	Color ambient_contribution = Color::scalarMultiplication(effective_color, m.m_ambient);			// Compute ambient contribution.
	real cosine_light_normal(Vector::dotProduct(lightv, this->normal));			// Cosine of angle between light vector and object surface normal. 
	if (cosine_light_normal < 0 || (shadow == true)) {														// If negative, the light is "in the other side" of the surface.
		diffuse = Color(0, 0, 0);
		specular = Color(0, 0, 0);
		//std::cout << "En SOMBRA!\n";
	}
	else {																				// If light and normal are in the same side...
		Color d = Color::scalarMultiplication(effective_color, (m.m_diffuse * cosine_light_normal));	// Compute diffusion contribution.		
		diffuse = Color(d.m_r, d.m_g, d.m_b);
		Vector lightv_n(Vector::negate(lightv)); 										// Compute the reflection of the light vector towards
		Vector reflectv(Vector::reflected(lightv_n, this->normal));							// the viewer ("eye" vector).
		real reflect_dot_eye(Vector::dotProduct(reflectv, this->eye));
		if (reflect_dot_eye <= 0) {				// If negative, there is no specular reflection contribution.
			specular = Color(0, 0, 0);
		}
		else {
			real factor(pow(reflect_dot_eye, m.m_shininess));										// If positive compute the specular contribution.
			Color s = Color::scalarMultiplication(light.m_intensity, (m.m_specular * factor));
			specular = Color(s.m_r, s.m_g, s.m_b);
		}
	}
	Color kodak = Color::addColors(ambient_contribution, (Color::addColors(diffuse, specular)));		// Add and return all reflection contributions.
	return kodak;
}

Color Tracer::shade_hit(World& w) {
		Color color;
		if (this->m_hit3D < 0 || this->m_hit3D == high_number) {
			color = black;
			return color;
		}
		else {
			std::unordered_map<real, Shape*>::const_iterator got = ix_shape_map.find(this->m_hit3D);
			Shape* s = got->second;
			bool shadow = shadowed(w, this->hit3D);
			color = lighting(this->s->b_material, w.m_light_source, shadow);
			return color;
		}
}

Color Tracer::color_at(World& w, Ray& r) {
		intersect_world(w, r);
		prepare_computations(r, w);
		if ( (this->m_hit3D == high_number)) { 
			return black; }
		else {
			this->ix_points.clear();				// Shadow rays use the same data structure as normal rays
			Color color = shade_hit(w);				// so the IXs vectors must be cleared.
			return color;
		}
}

void Tracer::render(Camera& camera, World& world) {
		std::ofstream tcfile;
		tcfile.open("test.ppm");
		tcfile << "P3\n" << camera.m_hsize << ' ' << camera.m_vsize << '\n' << 255 << '\n';
		Canvas image = Canvas(camera.m_hsize, camera.m_vsize);
		int vs = camera.m_vsize;
		int hs = camera.m_hsize;
		for (int y = 0; y < vs; y++) {
			for (int x = 0; x < hs; x++) {
				Ray ray = camera.ray_for_pixel(x, y);
				Color color = this->color_at(world, ray);
				this->ix_points.clear();				// if not clear the tracer array keeps adding records.
				this->ix_shape_map.clear();
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
	}
