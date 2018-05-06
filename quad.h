#ifndef _QUAD_H_
#define _QUAD_H_

#include "vec3.h"
#include "triangle.h"
#include "linear.h"

namespace Renzoku {
class Quad : public Shape {
public:
    Quad(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2, const Vec3 &_p3); 
    Quad(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2, const Vec3 &_p3,
         const Vec2 &uv0, const Vec2 &uv1, const Vec2 &uv2, const Vec2 &uv3); 
	
    inline bool hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const;
    inline bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;
    inline virtual void fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const;

    inline Float area() const;
    
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;  
	virtual Float pdf(const Vec3 &p) const;

    virtual void get_triangles(vector<Triangle> &list) const;
    inline virtual Triangle *get_triangle_pointer(int triangle_index) const;
    inline virtual int  get_triangle_count() const;
    inline virtual void get_triangle_coordinates(vector<Float> &coords) const;
    inline virtual void get_vertex_data(vector<GlVertex> &vertices) const;
    
    inline virtual Vec3 normal(const Vec3 &p) const;
    inline virtual Vec3 tangent() const;
    inline virtual Vec3 shading_normal(const Vec3 &p) const;
    inline virtual Vec2 texture_uv(const Vec3 &p) const;

    inline virtual void get_vertex_positions(vector<Vec3> &positions) const;    
    inline virtual void set_vertex_normal(int vertex_idx, const Vec3 &normal);
        
    inline virtual void set_surface_index(int idx);
    Quad& operator=(const Quad &q);

protected:
	inline virtual void compute_bounding_box();
    inline virtual void compute_bounding_sphere();

public:
    Vec3 p0, p1, p2, p3;
    Triangle t0, t1;
    
    Float quad_area;
    Float t0_ratio; 
    Triangle *first, *second;
};

inline Float Quad::area() const {
    return quad_area;
}

inline bool Quad::hit(const Ray &r, Float tmin, Float tmax, Float time) const {    
    if (first->hit(r, tmin, tmax, time))
        return true;
    return second->hit(r, tmin, tmax, time);
}

inline bool Quad::hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const {
    if (first->hit(r, tmin, tmax, time, record)) 
        return true;    
    return second->hit(r, tmin, tmax, time, record);
}

inline void Quad::fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const {
    record.copy_geometry_hit(gh);
    ((Triangle *)gh.shape)->fill_hit_record(r, gh, record);
    record.shape = (Shape *)this;        
    record.normal = first->normal(r.org() + gh.t * r.dir());
    record.tangent = first->tangent();
}

inline void Quad::compute_bounding_box() {
	box = BoundingBox(p0, p1, p2, p3);
}

inline void Quad::compute_bounding_sphere() {
    sphere = BoundingSphere(p0, p1, p2);
    sphere.merge(p3);
}

inline int  Quad::get_triangle_count() const {
    return 2;
}

inline Triangle *Quad::get_triangle_pointer(int triangle_index) const {
    if (triangle_index == 0) 
        return (Triangle *)&t0;
    else if (triangle_index == 1)
        return (Triangle *)&t1;
    else return NULL;
}

inline void Quad::get_triangles(vector<Triangle> &list) const {
    list.push_back(t0);
    list.push_back(t1);
}

inline void Quad::get_triangle_coordinates(vector<Float> &coords) const {
    t0.get_triangle_coordinates(coords);
    t1.get_triangle_coordinates(coords);    
}

inline void Quad::get_vertex_data(vector<GlVertex> &vertices) const {
    t0.get_vertex_data(vertices);
    t1.get_vertex_data(vertices);
}

inline Vec3 Quad::normal(const Vec3 &p) const {
    return first->normal(p);
}

inline Vec3 Quad::tangent() const {
    return first->tangent();
}

inline Vec3 Quad::shading_normal(const Vec3 &p) const {
    Float a, b, c;
    Linear::solve33(p0.x(), p1.x(), p2.x(),
                    p0.y(), p1.y(), p2.y(),
                    p0.z(), p1.z(), p2.z(),
                    p.x(), p.y(), p.z(), 
                    a, b, c);
    if (a >= 0.0f && b >= 0.0f && c >= 0.0f && a + b + c <= 1.0f)
        return unit_vector(a * t0.n0 + b * t0.n1 + c * t0.n2);

    return t1.shading_normal(p);
}

inline Vec2 Quad::texture_uv(const Vec3 &p) const {
    Float a, b, c;
    Linear::solve33(p0.x(), p1.x(), p2.x(),
                    p0.y(), p1.y(), p2.y(),
                    p0.z(), p1.z(), p2.z(),
                    p.x(), p.y(), p.z(), 
                    a, b, c);
    if (a >= 0.0f && b >= 0.0f && c >= 0.0f && a + b + c <= 1.0f)
        return a * t0.uv0 + b * t0.uv1 + c * t0.uv2;

    return t1.texture_uv(p);
}


inline void Quad::get_vertex_positions(vector<Vec3> &positions) const {
    positions.push_back(p0);
    positions.push_back(p1);
    positions.push_back(p2);
    positions.push_back(p3);
}

inline void Quad::set_surface_index(int idx) {
    this->surface_index = idx;
    t0.set_surface_index(idx);
    t1.set_surface_index(idx);
}

} // end namespace Renzoku

#endif 

