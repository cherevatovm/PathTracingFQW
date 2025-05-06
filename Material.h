#ifndef MATERIAL_H
#define MATERIAL_H
#include "Vec.h"

enum Refl_type { DIFF, SPEC, REFR, ROUGH };

struct Material {
	Refl_type brdf;
	Vec color;
	double refr_ind;
	double roughness;

	Material(const Vec& color_ = Vec(0.7, 0.7, 0.7), Refl_type brdf_ = DIFF,
		double refr_ind_ = 1, double roughness_ = 0) : 
		brdf(brdf_), color(color_), 
		refr_ind(refr_ind_), roughness(roughness_) {}
};

Material plastic = {
	Vec(0.7, 0.7, 0.7),
	DIFF
};
Material mirror = {
	Vec(0.85, 0.85, 0.85),
	SPEC
};
Material glass = {
	Vec(0.85, 0.85, 0.85),
	REFR,
	1.5
};
Material metal = {
	Vec(0.75, 0.75, 1),
	ROUGH,
	0.1,
	0.25
};

std::vector<Material> materials = { plastic, mirror, glass, metal };

#endif
