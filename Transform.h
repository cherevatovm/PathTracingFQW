#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "SquareMatrix.h"
#include "Vec.h"
#include "Bounds.h"

class Transform {
private:
	SquareMatrix<4> matr, inv_matr;

public:
	Transform() = default;

	Transform(const SquareMatrix<4>& m, const SquareMatrix<4>& inv_m)
		: matr(m), inv_matr(inv_m) {
	}

	bool operator==(const Transform& t) const { return t.matr == matr; }

	Vec apply_transform(const Vec& v, bool is_point = true) const {
		double x_new = matr.contents[0][0] * v.x + matr.contents[0][1] * v.y + matr.contents[0][2] * v.z + matr.contents[0][3];
		double y_new = matr.contents[1][0] * v.x + matr.contents[1][1] * v.y + matr.contents[1][2] * v.z + matr.contents[1][3];
		double z_new = matr.contents[2][0] * v.x + matr.contents[2][1] * v.y + matr.contents[2][2] * v.z + matr.contents[2][3];
		double w_new = matr.contents[3][0] * v.x + matr.contents[3][1] * v.y + matr.contents[3][2] * v.z + matr.contents[3][3];
		if (w_new == 1 || !is_point)
			return Vec(x_new, y_new, z_new);
		else
			return Vec(x_new / w_new, y_new / w_new, z_new / w_new);
	}

	Vec transform_normal(const Vec& n) const {
		return Vec(inv_matr.contents[0][0] * n.x + inv_matr.contents[1][0] * n.y + inv_matr.contents[2][0] * n.z,
			inv_matr.contents[0][1] * n.x + inv_matr.contents[1][1] * n.y + inv_matr.contents[2][1] * n.z,
			inv_matr.contents[0][2] * n.x + inv_matr.contents[1][2] * n.y + inv_matr.contents[2][2] * n.z);
	}

	Bounds3 transform_bounds(const Bounds3& aabb) const {
		Bounds3 aabb_t;
		for (int i = 0; i < 8; ++i)
			aabb_t = Bounds3::find_union(aabb_t, apply_transform(aabb.get_corner(i)));
		return aabb_t;
	}

	const SquareMatrix<4>& get_matr() const { return matr; }
	const SquareMatrix<4>& get_inverse_matr() const { return inv_matr; }
};

#endif