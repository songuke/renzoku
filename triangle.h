#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "shape.h"
#include "vec3.h"
#include "material.h"
#include "linear.h"

namespace Renzoku {
/**
 * Front face convention: CCW (same as OpenGL).
 */ 
class Triangle : public Shape {
public: 
    Triangle();
    void init(const Vec3 &p0, const Vec3& p1, const Vec3 &p2);

    Triangle(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2);
    Triangle(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2, 
             const Vec2 &uv0, const Vec2 &uv1, const Vec2 &uv2);

    bool hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;    
    virtual void fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const;

    inline virtual Vec3 normal(const Vec3 &p) const;
    inline virtual Vec3 tangent() const;
    inline virtual Vec3 shading_normal(const Vec3 &p) const;
    inline virtual Vec2 texture_uv(const Vec3 &p) const;

    inline Float area() const;

    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;
    virtual Float pdf(const Vec3 &p) const;
    
    inline virtual int  get_triangle_count() const;
        
    virtual void get_triangles(vector<Triangle> &list) const;
    virtual void get_triangle_coordinates(vector<Float> &coords) const;
    virtual Triangle *get_triangle_pointer(int triangle_index) const;
    virtual void get_vertex_data(vector<GlVertex> &vertices) const;
        
    virtual void get_vertex_positions(vector<Vec3> &positions) const;
    virtual void set_vertex_normal(int vertex_idx, const Vec3 &normal);

    virtual void set_tangent(const Vec3 &tangent);


protected:
    inline virtual void compute_bounding_box();
    inline virtual void compute_bounding_sphere();

public:
    Vec3 p0, p1, p2;
    Vec3 n0, n1, n2;        // vertex normals
    Float tri_area;
    Vec3 n, u;              // store facet normal and tangent
    Vec2 uv0, uv1, uv2;     // texture mapping
}; 

inline Float Triangle::area() const {
    return tri_area;
}

inline Vec3 Triangle::tangent() const {
    return u;
}

inline Vec3 Triangle::normal(const Vec3 &p) const {
    return n;
}

inline Vec3 Triangle::shading_normal(const Vec3 &p) const {
    Float a, b, c;
    Linear::solve33(p0.x(), p1.x(), p2.x(),
                    p0.y(), p1.y(), p2.y(),
                    p0.z(), p1.z(), p2.z(),
                    p.x(), p.y(), p.z(), 
                    a, b, c);

    return unit_vector(a * n0 + b * n1 + c * n2);
}

inline Vec2 Triangle::texture_uv(const Vec3 &p) const {
    Float a, b, c;
    Linear::solve33(p0.x(), p1.x(), p2.x(),
                    p0.y(), p1.y(), p2.y(),
                    p0.z(), p1.z(), p2.z(),
                    p.x(), p.y(), p.z(), 
                    a, b, c);

    return a * uv0 + b * uv1 + c * uv2;
}

inline void Triangle::compute_bounding_box() {
	box = BoundingBox(p0, p1, p2);
}

inline void Triangle::compute_bounding_sphere() {
    sphere = BoundingSphere(p0, p1, p2);
}

inline int  Triangle::get_triangle_count() const {
    return 1;
}

inline Triangle *Triangle::get_triangle_pointer(int triangle_index) const {
    if (triangle_index == 0) 
        return (Triangle *)(this);
    return NULL;
}

inline void Triangle::get_triangles(vector<Triangle> &list) const {
    list.push_back(*this);
}

} // end namespace

#endif
