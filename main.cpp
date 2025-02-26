#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <ctime>
#include <chrono>
#include <random>
#include <vector>
#include "Vec.h"
#include "Ray.h"
#include "Shape.h"
#include "Rand_om.h"
#include "stb_image_write.h"

using namespace std::chrono;

std::vector<Sphere*> light_sources = { 
	new Sphere(5.5, Vec(50, 81.6 - 15.5, 90), Vec(1, 1, 1) * 40, Vec(), DIFF),
	//new Sphere(2.75, Vec(70, 81.6 - 30.5, 81.6), Vec(0, 1, 1) * 20,  Vec(), DIFF),
    //new Sphere(1.375, Vec(20, 81.6 - 20, 81.6), Vec(1, 1, 0) * 40,  Vec(), DIFF)
};

// radius, position, emission, color, material
//Sphere objects[] = {
std::vector<Shape*> objects = {
	new Triangle(Vec(0, 0, 170), Vec(0, 82.5, 170), Vec(), Vec(), Vec(0.75, 0.25, 0.25), DIFF), // Left wall
	new Triangle(Vec(0, 82.5, 170), Vec(0, 82.5, 0), Vec(), Vec(), Vec(0.75, 0.25, 0.25), DIFF), // Left wall
	
	new Triangle(Vec(99.5, 0, 0), Vec(99.5, 82.5, 0), Vec(99.5, 0, 170), Vec(), Vec(0.25, 0.25, 0.75), DIFF), // Right wall
	new Triangle(Vec(99.5, 82.5, 0), Vec(99.5, 82.5, 170), Vec(99.5, 0, 170), Vec(), Vec(0.25, 0.25, 0.75), DIFF), // Right wall
	
	new Triangle(Vec(), Vec(0, 82.5, 0), Vec(99.5, 0, 0), Vec(), Vec(0.25, 0.75, 0.25), DIFF), // Back wall
	new Triangle(Vec(0, 82.5, 0), Vec(99.5, 82.5, 0), Vec(99.5, 0, 0), Vec(), Vec(0.25, 0.75, 0.25), DIFF), // Back wall
	
	new Triangle(Vec(0, 0, 170), Vec(0, 82.5, 170), Vec(99.5, 0, 170), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Front wall
	new Triangle(Vec(0, 82.5, 170), Vec(99.5, 82.5, 170), Vec(99.5, 0, 170), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Front wall
	
	new Triangle(Vec(0, 0, 170), Vec(), Vec(99.5, 0, 170), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Floor
	new Triangle(Vec(), Vec(99.5, 0, 0), Vec(99.5, 0, 170), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Floor

	new Triangle(Vec(0, 82.5, 0), Vec(0, 82.5, 170), Vec(99.5, 82.5, 0), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Ceiling
	new Triangle(Vec(0, 82.5, 170), Vec(99.5, 82.5, 170), Vec(99.5, 82.5, 0), Vec(), Vec(0.75, 0.75, 0.75), DIFF), // Ceiling
	
	new Sphere(16.5, Vec(73, 16.5, 95), Vec(), Vec(1, 1, 1), REFR), // Glass sphere
	new Sphere(20.5, Vec(33, 20.5, 65), Vec(), Vec(1, 1, 1), DIFF), // Matte sphere 
	
	light_sources[0], //light_sources[1], light_sources[2]
};


void create_orthonorm_sys(const Vec& v1, Vec& v2, Vec& v3) {
	// Projection to y = 0 plane and normalized orthogonal vector construction
	if (std::abs(v1.x) > std::abs(v1.y)) {
		double inv_len = 1.0 / sqrt(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * inv_len, 0.0, v1.x * inv_len);
	}
	// Projection to x = 0 plane and normalized orthogonal vector construction
	else {
		double inv_len = 1.0 / sqrt(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0, v1.z * inv_len, -v1.y * inv_len);
	}
	v3 = v1 % v2;
}

Vec sample_hemisphere(double u1, double u2) {
	const double r = sqrt(1.0 - u1 * u1);
	const double phi = 2 * M_PI * u2;
	return Vec(cos(phi) * r, sin(phi) * r, u1);
}

inline double clamp01(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int to_int(double x) { return int(pow(clamp01(x), 1 / 2.2) * 255 + 0.5); }

inline bool intersect_scene(const Ray& r, double& t, int& id) {
	double d;
	double inf = t = LDBL_MAX;
	for (int i = 0; i < objects.size(); ++i) {
		if ((d = objects[i]->intersect(r)) && d < t) {
			t = d;
			id = i;
		}
	}
	return t < inf;
}

Vec path_tracing(const Ray& r, int depth, unsigned short* Xi, int E = 1) {
 	if (depth > 100) return Vec();

	double t; // Distance to intersection 
	int id = 0; // ID of intersected object 
	if (!intersect_scene(r, t, id)) return Vec(); // Return black if there is no intersection 

	const Shape* hit_obj = objects[id];
	Vec hit_point = r.orig + r.dir * t;
	Vec n = hit_obj->get_normal(hit_point);
	Vec nl = n.dot_prod(r.dir) < 0 ? n : n * -1;
	Vec color = hit_obj->c;
	double rr_prob = std::max(color.x, std::max(color.y, color.z));

	// Russian Roulette for path termination
	if (++depth > 5 || !rr_prob) {
		if (erand48(Xi) < rr_prob)
			color = color * (1 / rr_prob);
		else
			return hit_obj->e * E;
	}

	// Diffuse BRDF
	if (hit_obj->refl == DIFF) {
		double r1 = 2 * M_PI * erand48(Xi);
		double r2 = erand48(Xi);
		double r2s = sqrt(r2);

		Vec w = nl, u, v;
		create_orthonorm_sys(w, u, v);
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

		Vec e;
		int s1 = objects.size(), s2 = light_sources.size();
		for (int i = 0; i < s2; ++i) {
			// Create orthonormal coord system and sample direction by solid angle
			Vec sw = (light_sources[i]->p - hit_point).norm(), su, sv;
			create_orthonorm_sys(sw, su, sv);

			double cos_a_max = sqrt(1 - light_sources[i]->rad * light_sources[i]->rad / 
				(hit_point - light_sources[i]->p).dot_prod(hit_point - light_sources[i]->p));
			double eps1 = erand48(Xi), eps2 = erand48(Xi);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * M_PI * eps2;

			Vec samp_dir = (su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a).norm();

			bool is_clear = intersect_scene(Ray(hit_point, samp_dir), t, id);
			// shoot shadow rays
			if (is_clear && (id == i + s1 - s2)) {
				double omega = 2 * M_PI * (1 - cos_a_max);
				e = e + color.mult(light_sources[i]->e * samp_dir.dot_prod(nl) * omega) * (1 / M_PI); // 1/PI for BRDF
			}
		}

		return hit_obj->e * E + e + color.mult(path_tracing(Ray(hit_point, d), depth, Xi, 0));
	}
	// Specular BRDF
	else if (hit_obj->refl == SPEC)
		return hit_obj->e + color.mult(path_tracing(Ray(hit_point, r.dir - n * n.dot_prod(r.dir) * 2), depth, Xi)); // Angle of incidence == angle of reflection

	// Refractive BRDF
	Ray refl_ray(hit_point, r.dir - n * n.dot_prod(r.dir) * 2);
	bool into = n.dot_prod(nl) > 0; // Check if a ray goes in or out 
	double refr_ind1 = 1, refr_ind2 = 1.5;
	double refr_ratio = into ? refr_ind1 / refr_ind2 : refr_ind2 / refr_ind1;
	double cos_incid_angle = r.dir.dot_prod(nl);
	double cos2t = 1 - refr_ratio * refr_ratio * (1 - cos_incid_angle * cos_incid_angle);

	if (cos2t < 0) 
		return hit_obj->e + color.mult(path_tracing(refl_ray, depth, Xi)); // Total internal reflection 

	Vec tdir = (r.dir * refr_ratio - n * ((into ? 1 : -1) * (cos_incid_angle * refr_ratio + sqrt(cos2t)))).norm();
	double a = refr_ind2 - refr_ind1;
	double b = refr_ind2 + refr_ind1;
	double c = 1 - (into ? -cos_incid_angle : tdir.dot_prod(n));

	double F0 = a * a / (b * b);
	double Fr = F0 + (1 - F0) * pow(c, 5);
	double Tr = 1 - Fr;

	double prob = 0.25 + 0.5 * Fr;
	double refl_prob = Fr / prob, trans_prob = Tr / (1 - prob);

	Vec result;
	if (depth > 2) {
		// Russian roulette
		if (erand48(Xi) < prob)
			result = path_tracing(refl_ray, depth, Xi) * refl_prob;
		else
			result = path_tracing(Ray(hit_point, tdir), depth, Xi) * trans_prob;
	}
	else
		result = path_tracing(refl_ray, depth, Xi) * Fr + path_tracing(Ray(hit_point, tdir), depth, Xi) * Tr;

	return hit_obj->e + color.mult(result);
}

int main(int argc, char* argv[]) {
	setlocale(LC_ALL, "RUSSIAN");

	int width = 1024, height = 768, samples = 32 / 4;
	Vec* image = new Vec[width * height], result;
	
	Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); 
	Vec camera_x = Vec(width * 0.5135 / height),
	camera_y = (camera_x % camera.dir).norm() * 0.5135;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic, 1) private(result)       // Use OpenMP 
	for (int y = 0; y < height; y++) {                       // Go through image rows 
		fprintf(stderr, "\rRendering (%d samples per pixel): %5.2f%%", samples * 4, 100.0 * y / (height - 1));
		for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < width; x++) // Go through image cols 
			// 2x2 subpixel rows 
			for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++)
				// 2x2 subpixel cols 
				for (int sx = 0; sx < 2; sx++, result = Vec()) {
					for (int s = 0; s < samples; s++) {
						double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						Vec d = camera_x * (((sx + 0.5 + dx) / 2 + x) / width - 0.5) +
							camera_y * (((sy + 0.5 + dy) / 2 + y) / height - 0.5) + camera.dir;
						result = result + path_tracing(Ray(camera.orig + d * 140, d.norm()), 0, Xi) * (1.0 / samples);
					}
					image[i] = image[i] + Vec(clamp01(result.x), clamp01(result.y), clamp01(result.z)) * 0.25;
				}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	printf("\nВремя работы: %d ms;\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());

	for (int i = 0; i < objects.size(); ++i)
		delete(objects[i]);

	// Save the result to png image
	auto* res_image = new unsigned char[width * height * 3];
	for (int i = 0; i < width * height; i++) {
		res_image[i * 3] = to_int(image[i].x);
		res_image[i * 3 + 1] = to_int(image[i].y);
		res_image[i * 3 + 2] = to_int(image[i].z);
	}

	stbi_write_png("result.png", width, height, 3, res_image, width * 3);
}