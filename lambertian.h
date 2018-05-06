#ifndef _LAMBERTIAN_H_
#define _LAMBERTIAN_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "onb.h"
#include "bsdf.h"
#include "rgb.h"

namespace Renzoku {

class Lambertian : public Bsdf, public IDiffuse {
public:
    Lambertian(Rgb kd) : kd(kd) {}
    Lambertian(Float _r, Float _g, Float _b) : kd(Rgb(_r, _g, _b)) {}      

public:
    Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);

    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    
    virtual Rgb get_diffuse_component()     const { return kd; }
    virtual Rgb get_representative_color()  const { return kd; }

    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::LAMBERTIAN; }
    
protected:
    Rgb kd;
};

} // end namespace Renzoku

#endif 
