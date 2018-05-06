#ifndef _SPOT_LIGHT_H_
#define _SPOT_LIGHT_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "rgb.h"
#include "light.h"
#include "onb.h"

namespace Renzoku {

/**
 * Spot light model: 
 * http://msdn.microsoft.com/en-us/library/windows/desktop/bb174697(v=vs.85).aspx
 * Falloff ratio: 
 * http://msdn.microsoft.com/en-us/library/windows/desktop/bb172279%28v=vs.85%29.aspx
 *
 * phi and theta represent the total angle of the outer and inner cone.
 *
 */
class SpotLight : public Light {
public:
    SpotLight(const Vec3 &_p, const Vec3 &_n, Float phi, Float theta, Float falloff, const Rgb &_intensity);
    
    virtual Light::Type get_light_type() const {
        return Light::SPOT_LIGHT;
    }

    bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    Float area() const;
    Rgb power() const;
        
    virtual void sample(Random &rd, Vec3 &p, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const;
    virtual void sample(Random &rd, Vec3 &p, Float &pdf) const;

    virtual Float pdf(const Vec3 &p) const;
    virtual Float pdf(const Vec3 &p, const Vec3 &wo) const;

    void query(const Vec3 &wo, Rgb &radiance) const;
    virtual Rgb first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler = NULL) const;
    virtual Rgb irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const;

    Vec3 org() const;
    Vec3 direction() const;

protected:    
    Vec3 pos, normal;
    Rgb intensity;
    Float falloff, phi, theta;
    Onb onb;
    
    Float cone_size;
    Float cos_half_phi, cos_half_theta;
    Float diff_cos;
};

} // end namespace Renzoku

#endif

