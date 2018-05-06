#ifndef _SHADING_GEOMETRY_H_
#define _SHADING_GEOMETRY_H_

#include "local_geometry.h"

namespace Renzoku {

/**
 * This class handles geometry and shading normal consistency. 
 * 
 * TODO: transform into local coordinate for eval and sampling.
 */
class ShadingGeometry {
public:
    ShadingGeometry(const Receiver &r);
    ShadingGeometry(const Ray &r, const HitRecord &rec);
    ShadingGeometry(const Vec3 &wo, Material *material, const LocalGeometry &dg);

    Rgb eval(const Vec3 &wi);
    Rgb sample(Random &rd, Vec3 &wi, Float &pdf);
    Float pdf(const Vec3 &wi);

    Float cosine(const Vec3 &wi);

protected:
    Vec3 wo;
    Material *material;
    LocalGeometry dg;
};

} // end namespace

#endif