#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include "Vec.h"
#include "Ray.h"
#include "Bounds.h"
#include "Transform.h"

enum Refl_type { DIFF, SPEC, REFR };

class Shape {
public:
	Vec e, c;
	Refl_type refl;
	const Transform* render_obj, *obj_render;

	Shape() = default;

	Shape(const Vec& e_, const Vec& c_, Refl_type refl_,
		const Transform* render_obj_ = nullptr, 
		const Transform* obj_render_ = nullptr) :
		e(e_), c(c_), refl(refl_), render_obj(render_obj_), obj_render(obj_render_) {};

	virtual double intersect(const Ray& r) const = 0;
	
	virtual Bounds3 get_bounds() const = 0;

	virtual Vec get_normal(const Vec& hit_point) const = 0;

	virtual Vec shape_sample(double u0, double u1) const = 0;
};

class Sphere : public Shape {
public:
	Vec p;
	double rad;

	Sphere() = default;

	Sphere(double rad_, Vec p_, 
		Vec e_, Vec c_, Refl_type refl_,
		const Transform* render_obj_ = nullptr, 
		const Transform* obj_render_ = nullptr) : Shape(e_, c_, refl_, render_obj_, obj_render_),
		rad(rad_), p(p_) {};
	
	/*
	double intersect(const Ray& r) const override { // returns distance, 0 if there is no hit 
		Vec op = p - r.orig; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t, eps = 1e-4,
			b = op.dot_prod(r.dir),
			det = b * b - op.dot_prod(op) + rad * rad;

		if (det < 0) return 0;
		else det = sqrt(det);

		return (t = b - det) > eps ? t :
			((t = b + det) > eps ? t : 0);
	}
	*/

	Sphere& operator=(const Sphere& sph) {
		p = sph.p;
		rad = sph.rad;
		e = sph.e; c = sph.c;
		refl = sph.refl;
		render_obj = sph.render_obj;
		obj_render = sph.obj_render;
		return *this;
	}

	// Returns distance or 0 if there is no hit 
	double intersect(const Ray& r) const override {
		Vec op = r.orig - p; // Solve t^2 * d.d + 2 * t * (o - p).d + (o - p).(o - p) - R^2 = 0 
		double eps = 1e-4;
		double b = 2 * (op).dot_prod(r.dir);
		double c = op.dot_prod(op) - rad * rad;
		double discr = b * b - 4 * c;

		if (discr < 0) return 0;
		else discr = sqrt(discr);

		double t0 = -b + discr;
		double t1 = -b - discr;

		return (t1 > eps) ? (0.5 * t1) :
			((t0 > eps) ? (0.5 * t0) : 0);
	}

	Vec get_normal(const Vec& hit_point) const override {
		return (hit_point - p).norm();
	}

	Bounds3 get_bounds() const override {
		return //render_obj->transform_bounds(Bounds3(
			Bounds3(Vec(p.x - rad, p.y - rad, p.z - rad), 
			Vec(p.x + rad, p.y + rad, p.z + rad));
	}

	Vec shape_sample(double u0, double u1) const override {
		return Vec();
	}
};

class Triangle : public Shape {
public:
	Vec vert1, vert2, vert3;
	Vec normal;

	Triangle() = default;

	Triangle(Vec v1_, Vec v2_, Vec v3_,
		Vec e_, Vec c_, Refl_type refl_,
		const Transform* render_obj_ = nullptr,
		const Transform* obj_render_ = nullptr) : 
		Shape(e_, c_, refl_, render_obj_, obj_render_),
		vert1(v1_), vert2(v2_), vert3(v3_),
		normal(((v2_ - v1_) % (v3_ - v1_)).norm()) {
	}

	Triangle& operator=(const Triangle& t) {
		vert1 = t.vert1;
		vert2 = t.vert2;
		vert3 = t.vert3;
		normal = t.normal;
		e = t.e; c = t.c;
		refl = t.refl;
		render_obj = t.render_obj;
		obj_render = t.obj_render;
		return *this;
	}

	double intersect(const Ray& r) const override {
		const double eps = 0.00001;

		double denominator = normal.dot_prod(r.dir);
		if (fabs(denominator) < eps)
			return 0;

		Vec ao = r.orig - vert1;
		double t = normal.dot_prod(ao) / denominator;
		t = -t;
		if (t < eps)
			return 0;

		Vec intersect = r.orig + r.dir * t;

		Vec v0 = vert2 - vert1;
		Vec v1 = vert3 - vert1;
		Vec v2 = intersect - vert1;

		double d00 = v0.dot_prod(v0);
		double d01 = v0.dot_prod(v1);
		double d11 = v1.dot_prod(v1);
		double d20 = v2.dot_prod(v0);
		double d21 = v2.dot_prod(v1);

		double denom = d00 * d11 - d01 * d01;
		double v = (d11 * d20 - d01 * d21) / denom;
		double w = (d00 * d21 - d01 * d20) / denom;
		double u = 1.0 - v - w;

		if (u >= 0 && v >= 0 && w >= 0)
			return t;

		return 0;
	}

	Vec get_normal(const Vec& hit_point) const override { return normal; }

	Bounds3 get_bounds() const override { 
		return Bounds3::find_union(Bounds3(vert1, vert2), vert3); 
	}

	Vec shape_sample(double u0, double u1) const override {
		return Vec();
	}
};

#endif