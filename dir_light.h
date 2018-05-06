#ifndef _DIR_LIGHT_H_
#define _DIR_LIGHT_H_

#include "common.h"
#include "light.h"
#include "rgb.h"
#include "vec3.h"
#include "size.h"

namespace Renzoku {

/**
 * Directional light.
 */
class DirectionalLight : public Light {
public:
    DirectionalLight(const Vec3 &wi, const Rgb &emission);
    
    virtual Light::Type get_light_type() const {
        return Light::DIRECTIONAL_LIGHT;
    }

    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    virtual Float area() const;
    virtual Rgb power() const;
    
    inline Vec3 get_wi() const;
    virtual void query(Vec3 &wi, Rgb &radiance) const;

    /**
     * Use a disk to sample parallel rays (as described in PBRT). 
     * A plane can be used for easier stratified sampling, but its size is less optimal.
     */
    void  set_approximate_area(Scene *scene);

    virtual Rgb first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler = NULL) const;
    virtual Rgb irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const;

protected:
    Rgb color;
    Vec3 wi;                        // point from surface towards the universe

    Float approx_area;
};

inline Vec3 DirectionalLight::get_wi() const {
    return wi;
}

} // end namespace Renzoku

#endif
