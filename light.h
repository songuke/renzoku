#ifndef _LIGHT_H_
#define _LIGHT_H_

#include "common.h"
#include "named_object.h"
#include "base_object.h"
#include "vec3.h"
#include "local_geometry.h"

namespace Renzoku {

struct LightSamplingRecord {
    bool is_valid;
    Vec3 wo;

    Vec3 p, n;
    bool has_location;
    bool has_normal;

    LightSamplingRecord() : is_valid(true), has_location(false), has_normal(false) {}
};

class Light : public BaseObject, public NamedObject {
public:
    enum Type {
        POINT_LIGHT,
        SPOT_LIGHT,
        AREA_LIGHT,
        DIRECTIONAL_LIGHT,
        ENV_LIGHT,
        BRDF_POINT_LIGHT,
        PHOTON_LIGHT
    };

public:
    Light() : prob(0.f) {}
    virtual ~Light() {}
    
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const = 0;
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const = 0;

    virtual Float area() const = 0;
    virtual Rgb power() const = 0;
    
    /**
     * Return the probability that this light is selected.
     */
    inline virtual Float selection_pdf() const;
    inline virtual void set_selection_pdf(Float prob);
        
    virtual Light::Type get_light_type() const = 0;

    /**
     * Return the irradiance at a receiver by integrating over all points on the light.
     *
     * E(x) = integrate over hemisphere L(x, w) cos(w) dw
     *      = integrate over all points on light L(x, y) G(x, y) dAy
     */
    virtual Rgb irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const = 0;
        
    /**
     * Return the reflected radiance at the receiver due to direct light.
     *
     * Visibility must be checked.
     *
     * Directional sampler can be used to pass a custom sampler to select a new direction.
     * For example, it can be a sampler that is built based on photon map distribution.
     */
    virtual Rgb first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler = NULL) const = 0;

protected:
    Float prob;                 // probability that this light is chosen due to power sampling
};

inline Float Light::selection_pdf() const {
    return prob;
}

inline void Light::set_selection_pdf(Float p) {
    prob = p;
}

} // end namespace Renzoku

#endif

