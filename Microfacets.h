#ifndef MICROFACETS_H
#define MICROFACETS_H
#define _USE_MATH_DEFINES
#include "Vec.h"

Vec random_hemisphere_dir(unsigned short* Xi, const Vec& n, double roughness) {
    double r1 = 2.0 * M_PI * erand48(Xi);
    double r2 = erand48(Xi);
    double alpha = roughness * roughness;

    //double theta = atan(sqrt(-alpha * log(1.0 - r2)));
    double theta = atan(alpha * sqrt(r2) / sqrt(1.0 - r2));
    double phi = r1;

    Vec localDir = Vec(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    );

    // Transform to world coords
    Vec w = n, u, v;
    create_orthonorm_sys(w, u, v);
    return (u * localDir.x + v * localDir.y + w * localDir.z).norm();
}

double beckmann_distribution(const Vec& n, const Vec& h, double roughness) {
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double alpha = roughness * roughness;
    double alpha2 = alpha * alpha;
    double n_dot_h_sqr = n_dot_h * n_dot_h;
    double tan2 = (1.0 - n_dot_h_sqr) / n_dot_h_sqr;
    return std::exp(-tan2 / alpha2) / (M_PI * alpha2 * n_dot_h_sqr * n_dot_h_sqr);
}

double ggx_distribution(const Vec& n, const Vec& h, double roughness) {
    double alpha = roughness * roughness;
    double alpha2 = alpha * alpha;
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double denom = n_dot_h * n_dot_h * (alpha2 - 1.0) + 1.0;
    return alpha2 / (M_PI * denom * denom);
}

double geom_cook_torrance(const Vec& n, const Vec& v, const Vec& l, const Vec& h, double roughness) {
    double n_dot_v = std::max(n.dot_prod(v), 0.0);
    double n_dot_l = std::max(n.dot_prod(l), 0.0);
    double n_dot_h = std::max(n.dot_prod(h), 0.0);
    double v_dot_h = std::max(v.dot_prod(h), 0.0);
    double G1 = (2.0 * n_dot_h * n_dot_v) / v_dot_h;
    double G2 = (2.0 * n_dot_h * n_dot_l) / v_dot_h;
    return std::min(1.0, std::min(G1, G2));
}

double fresnel_schlick(double cosTheta, double F0) { return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);  }

double geom_schlick_ggx(double n_dot_v, double roughness) {
    double k = (roughness + 1.0) * (roughness + 1.0) / 8.0;
    return n_dot_v / (n_dot_v * (1.0 - k) + k + 1e-6); // +1e-6 to prevent divison by 0
}

double geom_smith(const Vec& n, const Vec& v, const Vec& l, double roughness) {
    double n_dot_v = std::max(n.dot_prod(v), 0.0);
    double n_dot_l = std::max(n.dot_prod(l), 0.0);
    double ggx1 = geom_schlick_ggx(n_dot_v, roughness);
    double ggx2 = geom_schlick_ggx(n_dot_l, roughness);
    return ggx1 * ggx2;
}

#endif
