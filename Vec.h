#ifndef VEC_H
#define VEC_H
#include <math.h>
#include "Utility.h"

constexpr float eps = std::numeric_limits<float>::epsilon();

struct Vec {
	double x, y, z;

	Vec(double x_ = 0, double y_ = 0, double z_ = 0): x(x_), y(y_), z(z_) {}
	void print() const { printf("%f, %f, %f\n", x, y, z); }
	Vec operator+(const Vec& v) const { return Vec(x + v.x, y + v.y, z + v.z); }
	Vec& operator+=(const Vec& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vec operator-(const Vec& v) const { return Vec(x - v.x, y - v.y, z - v.z); }
	Vec operator-() const { return Vec(-x, -y, -z); }
	Vec operator*(double scalar) const { return Vec(x * scalar, y * scalar, z * scalar); }
	Vec& operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}
	Vec operator%(const Vec& v) const {
		return Vec(
			diff_of_products(y, v.z, z, v.y),
			diff_of_products(z, v.x, x, v.z),
			diff_of_products(x, v.y, y, v.x));
	}

	double operator[](int i) const { return (i == 0) ? x : (i == 1) ? y : z; }
	double dot_prod(const Vec& v) const { return x * v.x + y * v.y + z * v.z; }
	Vec mult(const Vec& v) const { return Vec(x * v.x, y * v.y, z * v.z); }
	Vec permute(const Vec& new_pos) const { return Vec((*this)[(int)new_pos[0]], (*this)[(int)new_pos[1]], (*this)[(int)new_pos[2]]); }
	double length_sqrd() const { return x * x + y * y + z * z; }
	double length() const { return sqrt(length_sqrd()); }

	Vec& norm() {
		double coef = 1 / length();
		x *= coef;
		y *= coef;
		z *= coef;
		return *this;
	}

	int min_comp_ind() {
		double min_val = std::min(x, std::min(y, z));
		if (std::abs(x - min_val) < eps)
			return 0;
		else if (std::abs(y - min_val) < eps)
			return 1;
		else
			return 2;
	}

	int max_comp_ind() {
		double max_val = std::max(x, std::max(y, z));
		if (std::abs(x - max_val) < eps)
			return 0;
		else if (std::abs(y - max_val) < eps)
			return 1;
		else
			return 2;
	}

	inline static Vec min(const Vec& v1, const Vec& v2) {
		return Vec(std::min(v1.x, v2.x),
			std::min(v1.y, v2.y),
			std::min(v1.z, v2.z));
	}

	inline static Vec max(const Vec& v1, const Vec& v2) {
		return Vec(std::max(v1.x, v2.x),
			std::max(v1.y, v2.y),
			std::max(v1.z, v2.z));
	}

	inline static Vec abs(const Vec& v) { return Vec(std::abs(v.x), std::abs(v.y), std::abs(v.z)); }

	void rotate_y(double alpha, const Vec& center) {
		double xc = x - center.x;
		double zc = z - center.z;
		double nx = xc * cos(alpha) - zc * sin(alpha);
		x = nx + center.x;
		double nz = xc * sin(alpha) + zc * cos(alpha);
		z = nz + center.z;
	}
};

#endif
