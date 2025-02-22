#ifndef RAY_H
#define RAY_H
#include "Vec.h"

struct Ray {
	Vec orig, dir;
	Ray(Vec o_, Vec d_) : orig(o_), dir(d_) {}
};

#endif