#ifndef _WARD_H_
#define _WARD_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "onb.h"
#include "bsdf.h"
#include "rgb.h"

namespace Renzoku {

/* 
 * Reference: 
 *   Ward, Measuring and Modeling Anisotropic Reflection, SIGGRAPH 1992.
 *   Walter, Notes on the Ward BRDF, 2005. 
 */
class Ward : public Bsdf, public IDiffuse, public ISpecular {
public:
    Ward(Rgb _kd, Rgb _ks, Float _ax, Float _ay);    

public:
    Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);
    
    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    
    virtual Rgb     get_diffuse_component()     const { return kd; }
    virtual Rgb     get_specular_component()    const { return ks; }    
    virtual Vec2    get_glossy_coefficients()   const { return Vec2(ax, ay); }
    virtual Rgb     get_representative_color()  const {
        if (kd.luminance() > 0)
            return kd;
        else 
            return ks;
    }    

    virtual Float   get_gloss_exponential() const {
        return std::max(1.0f / ax, 1.0f / ay);
    }

    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::WARD; }

protected:    
    Rgb kd, ks;
    Float ax, ay;       // elliptical Gaussian parameters to control anisotropy

    Float rho_d, rho_s; // for sampling

    Float ay_ax;
    Float ax2, ay2;     // precomputed values 
    Float inv_four_pi_ax_ay;
};

} // end namespace Renzoku

#endif 
