#ifndef _LOCAL_GEOMETRY_H_
#define _LOCAL_GEOMETRY_H_

#include "math3.h"
#include "common.h"
#include "hitrecord.h"
#include "ray.h"

namespace Renzoku {
    
struct Receiver {
    Vec3 p;
    Vec3 n;              // for normal consistency check
    Vec3 wo;    
           
    Material *m;
    AreaLight *light;    

    Vec2 uv;

    Vec3 tangent;
    Vec3 shading_n;
    Onb shading_uvn;
    
    Receiver() {}

    Receiver(const Ray &r, const HitRecord &h) : 
        p(r.org() + h.t * r.dir()), n(h.normal), wo(-r.dir()), m(h.material), light(h.light), shading_n(h.shading_normal), tangent(h.tangent), uv(h.uv), 
        shading_uvn(h.uvn)
    {

    }
};

/**
 * Represent local surface for BRDF sampling.
 */
struct LocalGeometry {
    Vec3 p;       
    Vec3 n;                     // shading normal    
    Vec3 ng;                    // geometric normal
    Vec3 tangent;    
    Vec2 uv;
    Onb uvn;

    LocalGeometry() {}

    LocalGeometry(const Vec3 &p, const Vec3 &n, const Vec3 &ng, const Onb &uvn, const Vec3 &tangent, const Vec2 uv)
        : p(p), n(n), ng(ng), uvn(uvn), tangent(tangent), uv(uv)
    {
    }

    LocalGeometry(const Ray &r, const HitRecord &h) 
        : p(r.org() + h.t * r.dir()), n(h.shading_normal), ng(h.normal), tangent(h.tangent), uv(h.uv), 
        uvn(h.uvn)
    {
    }

    LocalGeometry(const Receiver &r) 
        : p(r.p), n(r.shading_n), ng(r.n), tangent(r.tangent), uv(r.uv), uvn(r.shading_uvn)
    {
    }
};


/**
 * A general interface to generate a direction.
 * 
 * This is useful when a receiver needs a direction, e.g., to evaluate direct illumination.
 * To sample a new direction, it is possible to use BRDF sampling, 
 * or more advance techniques such as photon map or VPL map sampling.
 */
struct DirectionSampler {
    virtual Rgb sample(Random &rd, const Receiver &r, Vec3 &wi, Float &pdf) = 0;
    virtual Float pdf(const Receiver &r, const Vec3 &wi) = 0;

    virtual bool is_valid() = 0;
    virtual bool is_singular() = 0;
};

struct Convert {
    
    static Receiver to_receiver(const LocalGeometry &dg, Material *m, const Vec3 &wo)
    {
        Receiver r;
        r.p = dg.p;
        r.n = dg.ng;
        r.wo = wo;
        r.m = m;
        r.light = NULL;
        r.shading_n = dg.n;
        r.tangent = dg.tangent;
        r.uv = dg.uv;
        r.shading_uvn = dg.uvn;
        return r;
    }

};

} // end namespace

#endif