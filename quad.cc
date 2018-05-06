#include "quad.h"

namespace Renzoku {
   
Quad::Quad(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2, const Vec3 &_p3) 
    : p0(_p0), p1(_p1), p2(_p2), p3(_p3), 
      t0(_p0, _p1, _p2), t1(_p0, _p2, _p3)
{
    quad_area = t0.area() + t1.area();
    t0_ratio = t0.area() / quad_area;

    // for test hit, we test the larger triangle first    
    if (t0_ratio >= 0.5f) {
        first = &t0;
        second = &t1;
    } else {
        first = &t1;
        second = &t0;
    }

    // tangent should be shared for both of the triangles
    Vec3 p01 = unit_vector(p0 - p1);
    t0.set_tangent(p01);
    t1.set_tangent(p01);

    compute_bounding_box();
}


Quad::Quad(const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2, const Vec3 &_p3,
           const Vec2 &uv0, const Vec2 &uv1, const Vec2 &uv2, const Vec2 &uv3)
	: p0(_p0), p1(_p1), p2(_p2), p3(_p3), 
      t0(_p0, _p1, _p2, uv0, uv1, uv2), 
      t1(_p0, _p2, _p3, uv0, uv2, uv3)      
{
    quad_area = t0.area() + t1.area();
    t0_ratio = t0.area() / quad_area;

    // for test hit, we test the larger triangle first    
    if (t0_ratio >= 0.5f) {
        first = &t0;
        second = &t1;
    } else {
        first = &t1;
        second = &t0;
    }

    // tangent should be shared for both of the triangles
    Vec3 p01 = unit_vector(p0 - p1);
    t0.set_tangent(p01);
    t1.set_tangent(p01);

    compute_bounding_box();
}

Quad& Quad::operator=(const Quad &q) {
    p0 = q.p0;
    p1 = q.p1;
    p2 = q.p2;
    p3 = q.p3;
    t0 = q.t0;
    t1 = q.t1;    
    quad_area = q.quad_area;
    t0_ratio = q.t0_ratio;  
    if (t0_ratio >= 0.5f) {
        first = &t0;
        second = &t1;
    } else {
        first = &t1;
        second = &t0;
    }
    return *this;
}

void Quad::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    // sample depending on each triangle area per total area
    Float t = rd(); 
    if (t <= t0_ratio) {
        Float pdf_t0;
        t0.sample(rd, p, n, pdf_t0);
        pdf = t0_ratio * pdf_t0;
    } else {
        Float pdf_t1;
        t1.sample(rd, p, n, pdf_t1);
        pdf = (1.0f - t0_ratio) * pdf_t1;
    }
}

Float Quad::pdf(const Vec3 &p) const {
    return t0_ratio * t0.pdf(p) + (1.0f - t0_ratio) * t1.pdf(p);
}

void Quad::set_vertex_normal(int vertex_idx, const Vec3 &normal) {
    switch (vertex_idx) {
    case 0: 
        t0.set_vertex_normal(0, normal);
        t1.set_vertex_normal(0, normal);
        break;
    case 1:
        t0.set_vertex_normal(1, normal);
        break;
    case 2:
        t0.set_vertex_normal(2, normal);
        t1.set_vertex_normal(1, normal);
        break;
    case 3:
        t1.set_vertex_normal(2, normal);
        break;
    default:
        break;
    }
}

} // end namespace Renzoku

