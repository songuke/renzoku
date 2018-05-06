#include "sphere.h"
#include "scene.h"

namespace Renzoku {
Sphere::Sphere(const Vec3 &_center, Float _rad)
    : center(_center), rad(_rad) {
    
    sphere_area = 4 * A_PI * rad * rad;

    compute_bounding_box();
}

bool Sphere::hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const {
    record.hit = false;
    Vec3 tmp = r.org() - center;    

    Float a = dot(r.dir(), r.dir());
    Float b = dot(r.dir(), tmp); // simplify *2 
    Float c = dot(tmp, tmp) - rad * rad;
    
    Float discriminant = b*b - a*c;
    if (discriminant >= 0) {
        discriminant = sqrt(discriminant);
        
        Float t = (-b - discriminant) / a;
        if (t < tmin)
            t = (-b + discriminant) / a;
        
        if (t < tmin || t > tmax)
            return false;

        record.t = t;
        record.hit = true;
        record.shape = (Shape *)this;
        return true;
    }
    return false;
}

bool Sphere::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    Vec3 tmp = r.org() - center;    

    Float a = dot(r.dir(), r.dir());
    Float b = dot(r.dir(), tmp); // simplify *2 
    Float c = dot(tmp, tmp) - rad * rad;
    
    Float discriminant = b*b - a*c;
    if (discriminant >= 0) {
        discriminant = sqrt(discriminant);
        
        Float t = (-b - discriminant) / a;
        if (t < tmin)
            t = (-b + discriminant) / a;

        return (t >= tmin && t <= tmax);
    }
    return false;
}

void Sphere::fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const {
    record.copy_geometry_hit(gh);

    Vec3 p = r.org() + record.t * r.dir();
    record.normal = normal(p);
    record.tangent = tangent();        
    record.shading_normal = record.normal;
}

Float Sphere::area() const {
    return sphere_area;
}

/**
 * Uniform sampling is good for both surface in and out of the sphere.
 */
void Sphere::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    // uniform sampling over the surface area of the sphere: p(theta, phi) = sin(theta) / 4 * PI
    // this sampling is very bad as most of the points on the sphere can be occluded by the sphere itself when we connect
    // the surface point to the sampled point. 
    Vec2 uv;
    uv.random(rd);
    Float phi = 2 * M_PI * uv[0];
    Float cos_theta = 1.0 - 2 * uv[1];
    Float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    
    n = Vec3(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta);
    p = center + rad * n;

    pdf = 1.0f / sphere_area;
}

Float Sphere::pdf(const Vec3 &p) const {
    Float dist = (p - center).length();
    if (dist < rad - ZERO_EPSILON || dist > rad + ZERO_EPSILON) {
        return 0.0f;        
    }
    return 1.0f / sphere_area;
}

/**
 * Importance sampling is only good for surface outside the sphere. 
 */
void Sphere::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const {
    // determine if the point is out of the sphere
    Vec3 x = patch.p;    
    if ((x - center).length() < rad) { // in the sphere, use uniform sampling
        sample(rd, p, n, pdf);
        return;
    }

    // importance sample the sphere with additional information from the receiver
    Vec2 d;
    d.random(rd);
        
    Float cos_theta_max = sqrt(1.0f - (rad * rad) / (x - center).squared_length());
    
    Float cos_theta = 1.0f - d[0] * (1.0f - cos_theta_max);
    Float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    Float phi = TWO_PI * d[1];
    
    Onb basis;
    Vec3 nn = unit_vector(center - x);    
    basis.init_from_n(nn);

    Vec3 cone_dir = basis.local_to_world(Vec3(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta));

    GeometryHit rec;
    Float tick = 0.f;
    // local intersection test - only test the ray with this shape
    if (! this->hit(Ray(x, cone_dir), RAY_TMIN, RAY_TMAX, tick, rec)) {        
        // numerical inaccuracy near to the grazing angle, return the point in the center direction
        p = center - cone_dir * rad;
        n = normal(p);
    } else {        
        p = x + cone_dir * rec.t;
        n = normal(p);

        if (dot(x - p, n) <= ZERO_EPSILON) { // very near to grazing angle
            // to avoid zero pdf
            p = center - cone_dir * rad;
            n = normal(p);
        }
    }

    // probablity of the area
    Float pdf_w = 1.0f / (TWO_PI * (1.0f - cos_theta_max));
    pdf = dot(unit_vector(x - p), n) / (x - p).squared_length() * pdf_w;     

    /*
    if (pdf == 0.0) {
        cout << "Zero pdf" << endl;
        Vec3 xp = x - p;
        Float dot2 = dot(xp, n);
        cout << "Dot2: " << dot2 << endl;
    }*/
}

Float Sphere::pdf(const Vec3 &p, const Receiver &patch) const {
    Float dist = (p - center).length();
    if (dist < rad - ZERO_EPSILON || dist > rad + ZERO_EPSILON) {
        return 0.0f;        
    }

    Vec3 x = patch.p;    
    if ((x - center).length() < rad) { // in the sphere, use uniform sampling
        return pdf(p);
    }

    Vec3 n = normal(p);

    Float cos_theta_max = sqrt(1.0f - (rad * rad) / (x - center).squared_length());

    Float pdf_w = 1.0f / (TWO_PI * (1.0f - cos_theta_max));
    return dot(unit_vector(x - p), n) / (x - p).squared_length() * pdf_w;
}

int Sphere::get_triangle_count() const {
    return 0;
}

void Sphere::get_triangle_coordinates(vector<Float> &coords) const {
    // NOTE: can provide discretization of the sphere    
}

void Sphere::get_vertex_data(vector<GlVertex> &vertices) const {

}

Vec3 Sphere::normal(const Vec3 &p) const {
    return unit_vector(p - center);
}

Vec3 Sphere::tangent() const {
    /*
    Vec3 i(1.0, 0.0, 0.0);
    Vec3 j(0.0, 1.0, 0.0);
    Vec3 n = normal(p);
    Vec3 t = cross(n, i);
    if (t.length() < ONB_EPSILON)
        t = cross(n, j);
    return unit_vector(t);
    */
    return Vec3(1.0f, 0.0f, 0.0f);
}

}; // end namespace Renzoku
