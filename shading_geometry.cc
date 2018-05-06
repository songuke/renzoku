#include "shading_geometry.h"
#include "material.h"

namespace Renzoku {

ShadingGeometry::ShadingGeometry(const Receiver &r) : material(r.m), wo(r.wo), dg(LocalGeometry(r)) {
}

ShadingGeometry::ShadingGeometry(const Ray &r, const HitRecord &rec) : material(rec.material), wo(-r.dir()), dg(LocalGeometry(r, rec)) {
}

ShadingGeometry::ShadingGeometry(const Vec3 &wo, Material *material, const LocalGeometry &dg) : wo(wo), material(material), dg(dg) {
}

Float ShadingGeometry::cosine(const Vec3 &wi) {
    return fabs(dot(wi, dg.n));
}

Rgb ShadingGeometry::eval(const Vec3 &wi) {
    Vec3 ng = dg.ng;
    Vec3 ns = dg.n;

    Rgb brdf = DefaultRgb::black;
    
    bool strict_normals = false;
    // FIXME: this check makes the BRDF not symmetric anymore, which is bad! Disable for now.

    // geometry normal consistency check is performed here. 
    // we don't want to pass geometric normal into BRDF.
    //
    //                check by:
    //
    // specular:      cos(wi, ns) > 0 (BRDF)                      cos(wi, ng) > 0 (here)
    //                cos(wo, ns) > 0 (surface intersection)      cos(wo, ng) > 0 (here)
    // transmissive:  cos(wi, ns) < 0 (BRDF)                      cos(wi, ng) < 0 (here)
    //                cos(wo, ns) > 0 (surface intersection)      cos(wo, ng) > 0 (here)

    if (strict_normals) {
        if (dot(wo, ng) <= ZERO_EPSILON) {
            return DefaultRgb::black;
        }
        
        if (material->get_bsdf()->is_transmissive()) {
            if (dot(wi, ng) >= ZERO_EPSILON) {
                return DefaultRgb::black;
            }
        } else { // reflective
            if (dot(wi, ng) <= ZERO_EPSILON) {
                return DefaultRgb::black;
            }
        }
    } 
    brdf = material->eval(dg, wo, wi);    
    return brdf;
}

Rgb ShadingGeometry::sample(Random &rd, Vec3 &wi, Float &pdf) {    
    return material->sample(rd, dg, wo, wi, pdf);
}

Float ShadingGeometry::pdf(const Vec3 &wi) {
    return material->pdf(dg, wo, wi);
}

} // end namespace