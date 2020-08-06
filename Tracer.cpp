#include "Tracer.h"
#include <iostream>
#pragma omp parallel

using namespace rt;

Tracer::Tracer() {}

void Tracer::intersect_world(const World& w,const  Ray& ray) {
	Ray ray2(ray.m_origin, ray.m_direction);					// Ray must be transformed by inverse of shape transform matrix.
	for (int i = 0; i < w.m_objects.size(); i++) {
		Ray r = Ray::transform(ray2, w.m_objects[i]->b_inv_tx);
		w.m_objects[i]->intersection(r, this->ix_points, this->ix_shape_map);
		// refraction support map in test chapter 11
		// maps intersected object with their ID for easy retrieval.
		this->shape_id_map.insert({ w.m_objects[i]->b_id, w.m_objects[i] });
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
	sort_ix(this->ix_points);							// sort hit points.
	this->m_hit3D = min_hit_point(this->ix_points);		// return minimum positive hitpoint
	if (this->m_hit3D == high_number) return;			// if min+ hitpoint is high number there is no hit!
	std::unordered_map<real, Shape*>::const_iterator mapper = ix_shape_map.find(this->m_hit3D);	// Stores object mapped to the hit location
	this->s = mapper->second;															// according to hash map data structure.
	Matrix transform(this->s->b_inv_tx);				// Get the transform matrix of hit object.
	Vector eye = Vector::negate(ray.m_direction);		// Assigns eye vector the negative of ray direction.
	Point hitpoint_3d(Ray::position(ray, this->m_hit3D));		// Computes 3D coordinate of hit point as hitPoint = ray.orig + ray.dir * hit_location.
	s->normal_at(hitpoint_3d,this->s->b_inv_tx,this->normal);	// Normal at hitpoint after shape is transformed.
	if (Vector::dotProduct(normal, eye) < 0) {			// If normal and eye vector dot product is negative
		ray.m_inside = true;								// then the "viewer" is inside the object.
		normal = Vector::negate(normal);
	}
	else {
		ray.m_inside = false;
	}
	Point over_point = Point::addPoints(hitpoint_3d,			// compensate for shadows computation.
		Point::scalarMultiplication
		(Point(normal.m_x, normal.m_y, normal.m_z), epsilon));
	hitpoint_3d = over_point;
	this->hit3D = hitpoint_3d;									// Set hit point, eye, and normal vectors
	this->eye = eye;
	this->normal = normal;
	this->m_reflectv = Vector::reflected(ray.m_direction, this->normal); // Reflection vector at hit point.


	// Refraction implementation stage chapter 11
	for (int i = 0; i < this->ix_points.size(); i++) {
		if (this->ix_points[0] != infinite) {					// The code below does not make sense if 
			int verga = 0;										// the hit point is in the infinite!!!
			std::unordered_map<real, Shape*>::const_iterator mapper = ix_shape_map.find(ix_points[i]);
			Shape* s = mapper->second;
			if (this->io.empty() == true) {
				this->io.push_back(s->b_id);
			}
			else {
				for (int j = 0; j < io.size(); j++) {
					if (io[j] == s->b_id) {
						verga = io[j];
						this->io.erase(io.begin() + j);
						break;
					}
				}
			}
			if (verga == s->b_id) {
				int last_ix = this->io.back();
				std::unordered_map<int, Shape*>::const_iterator mapita = shape_id_map.find(last_ix);
				Shape* last_shape = mapita->second;
				fill_rix(this->n_ix, *last_shape);
				continue;
			}
			this->io.push_back(s->b_id);
			int last_ix = this->io.back();
			std::unordered_map<int, Shape*>::const_iterator mapita = shape_id_map.find(last_ix);
			Shape* last_shape = mapita->second;
			fill_rix(this->n_ix, *last_shape);
		}
		else continue;
	} 

	this->under_point = Point(hitpoint_3d.m_x - (normal.m_x * epsilon), // Point where refracted ray originates.
		hitpoint_3d.m_y - (normal.m_y * epsilon),
			hitpoint_3d.m_z - (normal.m_z * epsilon));

}

bool Tracer::shadowed(World& w, Point& p) {
	Vector v = Point::subPoints(w.m_light_source.m_source, p);
	real distance = Vector::magnitude(v);
	Vector direction = Vector::normalize(v);
	Ray r = Ray(p, direction);
	Tracer st = Tracer();
	st.intersect_world(w, r);
	st.m_hit3D = st.min_hit_point(st.ix_points);
	if (st.m_hit3D < 0 || st.m_hit3D == high_number) {
		return false;
	}
	else {
		if (st.m_hit3D < distance) return true; else return false;
	}
	st.ix_points.clear();
	st.ix_shape_map.clear();
}

Color Tracer::lighting(Material& m, Light& light, bool shadow) {
	Color diffuse(0, 0, 0);																// Diffuse reflection component.
	Color specular(0, 0, 0);															// Specular reflection component.
	Color effective_color;
	if (this->s->b_material.m_pattern.m_pattern_exist == true) {						// Compute color based on stripes if defined.
		effective_color = this->s->b_material.m_pattern.stripe_at_object(this->hit3D,this->s->b_inv_tx);
		// aqui llega el problema			std::cout << "test\n";
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
	else {																	// If light and normal are in the same side...
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

Color Tracer::shade_hit(World& w, real& remaining) {
		Color surface;
		Color reflected;
		Color polaroid;
		if (this->m_hit3D < 0 || this->m_hit3D == high_number) {
			surface = black;
			return surface;
		}
		else {
			std::unordered_map<real, Shape*>::const_iterator got = ix_shape_map.find(this->m_hit3D);
			Shape* s = got->second;
			bool shadow = shadowed(w, this->hit3D);
			surface = lighting(this->s->b_material, w.m_light_source, shadow);			
			// aqui llega el problema, lighting no retorna			std::cout << "test\n";
			if (this->s->b_material.m_reflective == 0) return surface;
			reflected = reflected_color( w, remaining);
			polaroid = Color::addColors(surface, reflected);
			return polaroid;
		}
}

Color Tracer::color_at(World& w, Ray& r, real& remaining) {
	this->intersect_world(w, r);
	this->prepare_computations(r, w);
		if ( (this->m_hit3D == high_number)) { 
			return black; }
		else {
			this->ix_points.clear();				// Shadow rays use the same data structure as normal rays
			this->io.clear();
			this->snell_sin.clear();
			this->n_ix.clear(); 
			Color color = this->shade_hit(w,remaining);				// so the IXs vectors must be cleared.
			// aqui llega el problema....std::cout << "first shade\n";
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
				real remaining = 2;
				Ray ray = camera.ray_for_pixel(x, y);
				Color color = this->color_at(world, ray, remaining);
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

Color Tracer::reflected_color( World& world, real& remaining) {
	real rr = remaining - 1;
	if (this->s->b_material.m_reflective == 0 || rr <= 0 || this->s->b_material.m_ambient == 1) {
		return black;
	}
	Tracer rt = Tracer();
	rt.hit3D = this->hit3D;
	rt.m_reflectv = this->m_reflectv;
	Ray reflect_ray(rt.hit3D, rt.m_reflectv);
	Color technicolor = Color::scalarMultiplication(rt.color_at(world, reflect_ray, rr)
		,this->s->b_material.m_reflective);
	rt.ix_points.clear();
	rt.ix_shape_map.clear();
	return technicolor;
}

real Tracer::snell(const real& n1, const real& n2) {
	this->n_ratio = 0;
	this->cos_i = 0;
	this->sin2_t = 0;
	n_ratio = n1 / n2;
	cos_i = Vector::dotProduct(this->eye, this->normal);
	sin2_t = pow(n_ratio, 2) * (1 - pow(cos_i, 2));
	return sin2_t;
}


void Tracer::fill_rix(std::vector<Refr_i>& n_ix, Shape& s) {
	Refr_i ri;
	if (n_ix.empty() == true) {			// if there are no intersections the 1st n1 = 1.0
		ri.n1 = 1.0;
		ri.n2 = s.b_material.m_refractive_index;
		this->snell_sin.push_back(snell(ri.n1, ri.n2));
	}
	else {
		ri.n1 = n_ix.back().n2;
		ri.n2 = s.b_material.m_refractive_index;
		this->snell_sin.push_back(snell(ri.n1, ri.n2));
	}
	if (this->io.size() == 1) {			// if there is one intersect object left then last n2 = 1.0
		ri.n1 = s.b_material.m_refractive_index;
		ri.n2 = 1.0;
		this->snell_sin.push_back(snell(ri.n1, ri.n2));
	}
	n_ix.push_back(ri);
} 

Color Tracer::refracted_color(World& world, real& remaining, const int& index) {
	real rr = remaining - 1; 
	if (this->s->b_material.m_transparecy == 0 || rr <= 0 || this->snell_sin[index] >= 1) {
		return black;
	}
	Tracer refract_tracer = Tracer(); 
	refract_tracer.hit3D = this->hit3D;
	refract_tracer.s = this->s;
	refract_tracer.s->b_inv_tx = this->s->b_inv_tx;
	real n_ratio = this->n_ix[index].n1 / this->n_ix[index].n2;
	real cos_t = sqrt(1.0 - this->snell_sin[index]);
	Vector direction = Vector::subVectors(
	Vector::scalarMultiplication(this->normal, (n_ratio * (this->cos_i - cos_t))),
	Vector::scalarMultiplication(this->eye, this->n_ratio));
	Ray refract_ray = Ray(this->under_point, direction); 
	Color palette = Color::scalarMultiplication(refract_tracer.color_at(world, refract_ray, rr)
		, this->s->b_material.m_transparecy); 
	refract_tracer.ix_points.clear();
	//refract_tracer.ix_shape_map.clear();
	refract_tracer.n_ix.clear();
	refract_tracer.io.clear();
	refract_tracer.snell_sin.clear();
		return palette;
}