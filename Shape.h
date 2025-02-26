#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include "Vec.h"
#include "Ray.h"
#include "Bounds.h"
#include "Transform.h"

const double eps = 1e-5;
enum Refl_type { DIFF, SPEC, REFR };

class Shape {
public:
	Vec emis, color;
	Refl_type refl;
	const Transform* render_obj, *obj_render;

	Shape() = default;

	Shape(const Vec& emis_, const Vec& color_, Refl_type refl_,
		const Transform* render_obj_ = nullptr, 
		const Transform* obj_render_ = nullptr) :
		emis(emis_), color(color_), refl(refl_), render_obj(render_obj_), obj_render(obj_render_) {};

	virtual double intersect(const Ray& r) const = 0;
	
	virtual Bounds3 get_bounds() const = 0;

	virtual Vec get_normal(const Vec& hit_point) const = 0;

	virtual Vec shape_sample(double u0, double u1) const = 0;
};

class Sphere : public Shape {
public:
	Vec center;
	double radius;

	Sphere() = default;

	Sphere(double radius_, Vec center_, 
		Vec emis_, Vec color_, Refl_type refl_,
		const Transform* render_obj_ = nullptr, 
		const Transform* obj_render_ = nullptr) : Shape(emis_, color_, refl_, render_obj_, obj_render_),
		radius(radius_), center(center_) {};
	
	/*
	// Returns distance or 0 if there is no hit
	double intersect(const Ray& r) const override {
		Vec op = p - r.orig; // Solve t^2 * d.d + 2 * t * (o - p).d + (o - p).(o - p) - R^2 = 0 
		double t, eps = 1e-4,
			b = op.dot_prod(r.dir),
			det = b * b - op.dot_prod(op) + radius * radius;

		if (det < 0) return 0;
		else det = sqrt(det);

		return (t = b - det) > eps ? t :
			((t = b + det) > eps ? t : 0);
	}
	*/

	Sphere& operator=(const Sphere& sph) {
		center = sph.center;
		radius = sph.radius;
		emis = sph.emis; color = sph.color;
		refl = sph.refl;
		render_obj = sph.render_obj;
		obj_render = sph.obj_render;
		return *this;
	}

	// Returns distance or 0 if there is no hit 
	double intersect(const Ray& r) const override {
		Vec op = r.orig - center; // Solve t^2 * d.d + 2 * t * (o - p).d + (o - p).(o - p) - R^2 = 0
		double b = 2 * (op).dot_prod(r.dir);
		double c = op.dot_prod(op) - radius * radius;
		double discr = b * b - 4 * c;

		if (discr < 0) return 0;
		else discr = sqrt(discr);

		double t0 = -b + discr;
		double t1 = -b - discr;

		return (t1 > eps) ? (0.5 * t1) :
			((t0 > eps) ? (0.5 * t0) : 0);
	}

	Vec get_normal(const Vec& hit_point) const override {
		return (hit_point - center).norm();
	}

	Bounds3 get_bounds() const override {
		return //render_obj->transform_bounds(Bounds3(
			Bounds3(Vec(center.x - radius, center.y - radius, center.z - radius), 
			Vec(center.x + radius, center.y + radius, center.z + radius));
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

	Triangle(Vec vert1_, Vec vert2_, Vec vert3_,
		Vec emis_, Vec color_, Refl_type refl_,
		const Transform* render_obj_ = nullptr,
		const Transform* obj_render_ = nullptr) : 
		Shape(emis_, color_, refl_, render_obj_, obj_render_),
		vert1(vert1_), vert2(vert2_), vert3(vert3_),
		normal(((vert2_ - vert1_) % (vert3_ - vert1_)).norm()) {
	}

	Triangle& operator=(const Triangle& t) {
		vert1 = t.vert1;
		vert2 = t.vert2;
		vert3 = t.vert3;
		normal = t.normal;
		emis = t.emis; color = t.color;
		refl = t.refl;
		render_obj = t.render_obj;
		obj_render = t.obj_render;
		return *this;
	}

	double intersect(const Ray& r) const override {
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

class Hexahedron {
public:
	std::vector<Triangle*> faces;

	Hexahedron(
		Vec center_,
		double radius_,
		double alpha_,
		Vec emis_, Vec color_,
		Refl_type refl_) {
		Vec p1{ center_.x + radius_, center_.y + radius_, center_.z + radius_ };
		Vec p7{ center_.x - radius_, center_.y - radius_, center_.z - radius_ };
		Vec p2{ p1.x, p1.y, p7.z };
		Vec p3{ p7.x, p1.y, p7.z };
		Vec p4{ p7.x, p1.y, p1.z };
		Vec p5{ p1.x, p7.y, p1.z };
		Vec p6{ p1.x, p7.y, p7.z };
		Vec p8{ p7.x, p7.y, p1.z };

		p1.rotate_y(alpha_, center_);
		p2.rotate_y(alpha_, center_);
		p3.rotate_y(alpha_, center_);
		p4.rotate_y(alpha_, center_);
		p5.rotate_y(alpha_, center_);
		p6.rotate_y(alpha_, center_);
		p7.rotate_y(alpha_, center_);
		p8.rotate_y(alpha_, center_);

		faces.push_back(new Triangle(p1, p2, p3, emis_, color_, refl_));
		faces.push_back(new Triangle(p3, p4, p1, emis_, color_, refl_));
		faces.push_back(new Triangle(p8, p7, p6, emis_, color_, refl_));
		faces.push_back(new Triangle(p6, p5, p8, emis_, color_, refl_));
		faces.push_back(new Triangle(p2, p6, p7, emis_, color_, refl_));
		faces.push_back(new Triangle(p7, p3, p2, emis_, color_, refl_));
		faces.push_back(new Triangle(p1, p5, p6, emis_, color_, refl_));
		faces.push_back(new Triangle(p6, p2, p1, emis_, color_, refl_));
		faces.push_back(new Triangle(p8, p7, p3, emis_, color_, refl_));
		faces.push_back(new Triangle(p3, p4, p8, emis_, color_, refl_));
		faces.push_back(new Triangle(p1, p4, p8, emis_, color_, refl_));
		faces.push_back(new Triangle(p8, p5, p1, emis_, color_, refl_));
	}
};

#endif