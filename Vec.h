#ifndef VEC_H
#define VEC_H
#include <math.h> 

struct Vec {
	double x, y, z;

	Vec(double x_ = 0, double y_ = 0, double z_ = 0): x(x_), y(y_), z(z_) {}
	
	Vec operator+(const Vec& v) const { return Vec(x + v.x, y + v.y, z + v.z); }
	Vec operator-(const Vec& v) const { return Vec(x - v.x, y - v.y, z - v.z); }
	Vec operator*(double scalar) const { return Vec(x * scalar, y * scalar, z * scalar); }
	Vec operator%(const Vec& v) const { return Vec(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	double operator[](int i) const { return (i == 0) ? x : (i == 1) ? y : z; }
	double dot_prod(const Vec& v) const { return x * v.x + y * v.y + z * v.z; }
	Vec mult(const Vec& v) const { return Vec(x * v.x, y * v.y, z * v.z); }
	double length() const { return sqrt(x * x + y * y + z * z); }

	Vec& norm() {
		double coef = 1 / length();
		x *= coef;
		y *= coef;
		z *= coef;
		return *this;
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
};

inline double clamp(double val, double low, double high) {
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

#endif
