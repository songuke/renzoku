#ifndef _PHONG_H_
#define _PHONG_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "onb.h"
#include "bsdf.h"
#include "rgb.h"

#include <cmath>
using namespace std;

namespace Renzoku {

/**
 * Modified Phong model:
 * f(x, wi, wo) = kd / pi + 
 *                ks (n + 2) / 2pi * cos(alpha)^n
 * 
 * Reference: 
 *   Lafortune, Using Modified Phong Reflectance Model for Physically Based Rendering, 1994.
 */
class ModifiedPhong : public Bsdf, public IDiffuse, public ISpecular {
public:
    /**
     * To avoid NaN, the shininess is automatically rounded to nearest integer.
     */
    ModifiedPhong(Rgb _kd, Rgb _ks, Float _n);

public:
    Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);
    
    virtual Rgb     eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    virtual Float   pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
        
    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::PHONG; }

    virtual Rgb     get_diffuse_component()     const { return kd; }
    virtual Rgb     get_specular_component()    const { return ks; }
    virtual Vec2    get_glossy_coefficients()   const { return Renzoku::Vec2(n, 0); }
    virtual Rgb     get_representative_color()  const {
        if (kd.luminance() > 0)
            return kd;
        else 
            return ks;
    }

    virtual Float   get_gloss_exponential() const {
        return n;
    }

protected:
    Float rho_d, rho_s;
    Rgb kd, ks;
    Float n;        // specular exponential term
};

} // end namespace Renzoku

#endif 
