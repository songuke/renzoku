#ifndef _SPHERE_H_
#define _SPHERE_H_

#include "common.h"
#include "shape.h"
#include "vec3.h"
#include "rgb.h"

namespace Renzoku {
class Sphere : public Shape {
public:
    Sphere(const Vec3 &_center, Float _rad);
    
    bool hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;
    virtual void fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const;

    Float area() const;    
    
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;
    virtual Float pdf(const Vec3 &p) const;

    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const;
    virtual Float pdf(const Vec3 &p, const Receiver &patch) const;

    virtual int  get_triangle_count() const;
    virtual void get_triangle_coordinates(vector<Float> &coords) const;
    virtual void get_vertex_data(vector<GlVertex> &vertices) const;

    virtual Vec3 normal(const Vec3 &p) const;
    virtual Vec3 tangent() const;
    inline virtual Vec3 shading_normal(const Vec3 &p) const;
    inline virtual Vec2 texture_uv(const Vec3 &p) const;

    virtual void get_vertex_positions(vector<Vec3> &positions) const {}    
    virtual void set_vertex_normal(int vertex_idx, const Vec3 &normal) {}

protected:
	inline virtual void compute_bounding_box();
    inline virtual void compute_bounding_sphere();

public:
    Vec3 center;
    Float rad;
    Float sphere_area;
};

inline Vec3 Sphere::shading_normal(const Vec3 &p) const {
    return normal(p);
}


inline Vec2 Sphere::texture_uv(const Vec3 &p) const {
    return Vec2();
}

inline void Sphere::compute_bounding_box() {
    Float diag = rad * M_SQRT2;
	box = BoundingBox(center - diag, center + diag);
}

inline void Sphere::compute_bounding_sphere() {
    sphere = BoundingSphere(center, rad);
}

}

#endif
