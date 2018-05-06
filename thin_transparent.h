#ifndef _THIN_TRANSPARENT_H_
#define _THIN_TRANSPARENT_H_

#include "rgb.h"
#include "bsdf.h"

namespace Renzoku {

/**
 * A thin transparent surface allows light to pass through without scattering the incoming light. 
 *
 * In the absence of reflection of the incoming ray due to Fresnel effect, 
 * this material is equivalent to Glass with medium index set to 1.
 *
 * This is a singular material and cannot be sampled.
 */
class ThinTransparent : public Bsdf {
public:
    ThinTransparent(const Rgb &color);

    virtual Bsdf::Type get_bsdf_type() const { return Bsdf::THIN_TRANSPARENT; }

    virtual Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf);
    virtual Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf, Rgb &reflectivity);
    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi);

    virtual bool is_singular()     const { return true; }
    virtual bool is_transmissive() const { return true; }
    Rgb get_representative_color() const { return color; }

protected:
    Rgb color;
};

} // end namespace Renzoku

#endif