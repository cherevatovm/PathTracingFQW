#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include "Vec.h"
#include "Ray.h"
#include "Bounds.h"
#include "Transform.h"
#include "Utility.h"

using uint = unsigned int;
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

class Mesh {
private:
	Transform* last_transform = nullptr;

	class MeshTriangle : public Shape {
	public:
		const Mesh* mesh;
		uint vert0_ind, vert1_ind, vert2_ind;
		uint vert0_norm_ind, vert1_norm_ind, vert2_norm_ind;
		Vec surface_normal;	

		MeshTriangle() = default;

		MeshTriangle(const Mesh* mesh_, Vec color_, Refl_type refl_, 
			const Transform* render_obj_ = nullptr, 
			const Transform* obj_render_ = nullptr) : 
			Shape(Vec(), color_, refl_, render_obj_, obj_render_), mesh(mesh_) {}

		MeshTriangle(const Mesh* mesh_,
			uint v1_ind_, uint v2_ind_, uint v3_ind_,
			Vec emis_, Vec color_, Refl_type refl_,
			uint v1_n_ind_ = UINT_MAX, uint v2_n_ind_ = UINT_MAX, uint v3_n_ind_ = UINT_MAX,
			const Transform* render_obj_ = nullptr, const Transform* obj_render_ = nullptr) :
			Shape(emis_, color_, refl_, render_obj_, obj_render_), mesh(mesh_),
			vert0_ind(v1_ind_), vert1_ind(v2_ind_), vert2_ind(v3_ind_),
			vert0_norm_ind(v1_n_ind_), vert1_norm_ind(v2_n_ind_), vert2_norm_ind(v3_n_ind_)
		{
			surface_normal = ((mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind]) %
				(mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind])).norm();
		}

		MeshTriangle& operator=(const MeshTriangle& t) {
			mesh = t.mesh;
			vert0_ind = t.vert0_ind;
			vert1_ind = t.vert1_ind;
			vert2_ind = t.vert2_ind;
			vert0_norm_ind = t.vert0_norm_ind;
			vert1_norm_ind = t.vert1_norm_ind;
			vert2_norm_ind = t.vert2_norm_ind;
			emis = t.emis; color = t.color;
			refl = t.refl;
			render_obj = t.render_obj;
			obj_render = t.obj_render;
			return *this;
		}

		void calc_surface_normal() {
			surface_normal = ((mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind]) %
				(mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind])).norm();
		}

		void set_vert_ind(uint i, uint new_val) { 
			i == 0 ? vert0_ind = new_val : (i == 1 ? vert1_ind = new_val : vert2_ind = new_val);
		}
		
		void set_norm_ind(uint i, uint new_val) { 
			i == 0 ? vert0_norm_ind = new_val : (i == 1 ? vert1_norm_ind = new_val : vert2_norm_ind = new_val); 
		}

		double intersect(const Ray& r) const override {
			double denominator = surface_normal.dot_prod(r.dir);
			if (fabs(denominator) < eps)
				return 0;

			Vec ao = r.orig - mesh->vertices[vert0_ind];
			double t = surface_normal.dot_prod(ao) / denominator;
			t = -t;
			if (t < eps)
				return 0;

			Vec intersect = r.orig + r.dir * t;

			Vec v0 = mesh->vertices[vert1_ind] - mesh->vertices[vert0_ind];
			Vec v1 = mesh->vertices[vert2_ind] - mesh->vertices[vert0_ind];
			Vec v2 = intersect - mesh->vertices[vert0_ind];

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

		Vec get_normal(const Vec& hit_point) const override { return surface_normal; }

		Bounds3 get_bounds() const override {
			return Bounds3::find_union(Bounds3(mesh->vertices[vert0_ind], mesh->vertices[vert1_ind]),
				mesh->vertices[vert2_ind]);
		}

		Vec shape_sample(double u0, double u1) const override {
			return Vec();
		}
	};

	Vec centroid() const {
		double x_sum = 0, y_sum = 0, z_sum = 0;
		for (const auto& vertex : vertices) {
			x_sum += vertex.x;
			y_sum += vertex.y;
			z_sum += vertex.z;
		}
		size_t count = vertices.size();
		return { x_sum / count, y_sum / count, z_sum / count };
	}

	void affine_transformation(const Transform& tr) {
		for (Vec& vert : vertices)
			vert = tr.apply_transform(vert);
		for (Vec& norm : vert_normals)
			norm = tr.transform_normal(norm);
		for (MeshTriangle* face : faces)
			face->calc_surface_normal();
	}

	void calc_vertex_normals() {
		vert_normals.clear();
		for (uint i = 0; i < vertices.size(); ++i) {
			std::list<Vec> triangle_normals;
			for (MeshTriangle* face : faces) {
				if (face->vert0_ind == i || face->vert1_ind == i || face->vert2_ind == i)
					triangle_normals.push_back(face->surface_normal);
			}

			Vec vert_norm;
			for (const Vec& normal : triangle_normals)
				vert_norm += normal;
			vert_norm *= (1.0 / triangle_normals.size());
			vert_normals.push_back(vert_norm);
		}
		for (MeshTriangle* face : faces) {
			face->vert0_norm_ind = face->vert0_ind;
			face->vert1_norm_ind = face->vert1_ind;
			face->vert2_norm_ind = face->vert2_ind;
		}
	}

	void process_face(MeshTriangle* face, const std::vector<std::string>& face_data) {
		uint vertex_ind, norm_ind, tex_coord_ind;
		for (int i = 0; i < 3; ++i) {
			std::istringstream iss(face_data[i]);

			iss >> vertex_ind;
			face->set_vert_ind(i, --vertex_ind);
			char ch1 = iss.peek();
			if (ch1 == '/') {
				iss.ignore();
				ch1 = iss.peek();
				if (ch1 == '/') {
					iss.ignore();
					iss >> norm_ind;
					face->set_norm_ind(i, --norm_ind);
				}
				else if (isdigit(ch1)) {
					iss >> tex_coord_ind;
					//res_vert.tex_coords = vert_tex_coords[--tex_coord_ind];
					ch1 = iss.peek();
					if (ch1 == '/') {
						iss.ignore();
						iss >> norm_ind;
						face->set_norm_ind(i, --norm_ind);
					}
				}
			}
		}
		face->calc_surface_normal();
	}

public:
	std::vector<MeshTriangle*> faces;
	std::vector<Vec> vertices;
	std::vector<Vec> vert_normals;
	Vec color;
	Refl_type refl;

	Mesh(Vec color_, Refl_type refl_) : color(color_), refl(refl_) {};

	Mesh(const std::vector<Vec>& vertices_, const std::vector<Vec>& vert_normals_, 
		const std::vector<uint>& indices_, Vec color_, Refl_type refl_) :
		vertices(vertices_), vert_normals(vert_normals_), color(color_), refl(refl_) {
		bool has_normals = vert_normals_.size() > 0;
		uint inc = has_normals ? 6 : 3;
		for (int i = 0; i < indices_.size(); i += inc) {
			if (has_normals) {
				faces.push_back(new MeshTriangle(this, indices_[i], indices_[i + 2], indices_[i + 4],
					Vec(), color_, refl_, indices_[i + 1], indices_[i + 3], indices_[i + 5]));
			}
			else
				faces.push_back(new MeshTriangle(this, indices_[i], indices_[i + 1], indices_[i + 2], Vec(), color_, refl_));
		}
		if (!has_normals)
			calc_vertex_normals();
	};

	void rotate_by_center(double theta, int rotate_index) { 
		affine_transformation(Transform::rotate_around_point(centroid(), theta, rotate_index)); 
	}
	
	void scale_by_center(double x, double y, double z) {
		affine_transformation(Transform::scale_around_point(centroid(), x, y, z));
	}

	void translate(const Vec& delta) {
		affine_transformation(Transform::translate(delta));
	}

	void load_model(const std::string& file_name) {
		std::ifstream file(file_name);
		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << file_name << std::endl;
			return;
		}

		std::string line;
		while (std::getline(file, line)) {
			std::istringstream iss(line);
			std::string type;
			iss >> type;
			if (type.empty() || type[0] == '#') {
				continue;
			}
			if (type == "v") {
				float x, y, z;
				iss >> x >> y >> z;
				vertices.emplace_back(x, y, z);
			}
			else if (type == "vn") {
				float x, y, z;
				iss >> x >> y >> z;
				vert_normals.emplace_back(x, y, z);
			}
			else if (type == "f") {
				auto spl = split(line);
				for (int i = 3; i < spl.size(); ++i) {
					MeshTriangle* cur_triangle = new MeshTriangle(this, color, refl);
					faces.push_back(cur_triangle);
					process_face(cur_triangle, std::vector<std::string>{ spl[1], spl[i - 1], spl[i] });
				}
			}
		}
		file.close();
	}

	static Mesh* create_hexahedron(Vec center_, double radius_, 
		double alpha_, Vec color_, Refl_type refl_) {
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

		return new Mesh(std::vector<Vec>{ p1, p2, p3, p4, p5, p6, p7, p8 },
			std::vector<Vec>(),
			std::vector<uint>{ 0, 1, 2, 2, 3, 0,
			7, 6, 5, 5, 4, 7, 1, 5, 6, 6, 2,
			1, 0, 4, 5, 5, 1, 0, 7, 6, 2,
			2, 3, 7, 0, 3, 7, 7, 4, 0 }, color_, refl_);
	}
};

#endif