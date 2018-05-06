#ifndef _POINT_LIGHT_H_
#define _POINT_LIGHT_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "rgb.h"
#include "light.h"

namespace Renzoku {

class PointLight : public Light {
public:
    PointLight(const Vec3 &_p, Rgb _intensity);
    
    virtual Light::Type get_light_type() const {
        return Light::POINT_LIGHT;
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
    Rgb radiance() const;

protected:
    Vec3 p;
    Rgb intensity;
};

} // end namespace Renzoku

#endif

