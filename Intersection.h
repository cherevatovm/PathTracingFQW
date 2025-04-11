#ifndef INTERSECTION_H
#define INTERSECTION_H
#include "Shape.h"

class Shape;

class Intersection {
public:
	double t;
	const Shape* object = nullptr;

	Intersection() = default;
	Intersection(double t_, const Shape* object_ = nullptr) :
		t(t_), object(object_) {
	}

	bool operator<(const Intersection& inters) const { return t < inters.t; }
};

class TriangleIntersection : public Intersection {
public:
	Vec baryc_coords;

	TriangleIntersection() = default;
	TriangleIntersection(double t_, const Shape* object_, Vec baryc_coords_) :
		Intersection(t_, object_), baryc_coords(baryc_coords_) {
	}
};

#endif