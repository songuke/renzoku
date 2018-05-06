#ifndef _AREA_LIGHT_H_
#define _AREA_LIGHT_H_

#include "common.h"
#include "rgb.h"
#include "light.h"
#include "surface.h"

namespace Renzoku {
class AreaLight : public Light {
public:
    AreaLight(Surface *_surface, const Rgb &emission);

    virtual Light::Type get_light_type() const {
        return Light::AREA_LIGHT;
    }

    inline virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    inline virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    inline virtual Float area() const;
    inline virtual Rgb power() const;
    
    /**
     * Return a point on the light
     */
    inline virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;

    /**
     * Return a point and an outgoing direction from the light.
     */
    inline virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf) const;
    inline virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const;
    
    /** 
     * Return the probability of sampling a point on the light.
     */
    inline virtual Float pdf(const Vec3 &p) const;

    /**
     * Return the probability of sampling a point and a direction from the light.
     */
    inline virtual Float pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo) const;
    inline virtual void  pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo, Float &pdf_p, Float &pdf_wo) const;

    /**
     * Given a point on the light, return probability of a direction
     */
    inline virtual Float pdf_wo(const Vec3 &p, const Vec3 &wo) const;

    /**
     * Consider receiver location and sample.
     */
    inline virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const;
    inline virtual Float pdf(const Vec3 &p, const Receiver &patch) const;    

    /**
     * Query the out-going radiance at a point on the light. Useful for cases when some ray hits the light.
     */
    virtual void query(const Vec3 &p, const Vec3 &n, const Vec3 &wo, Rgb &radiance) const;
        
    inline Surface *get_surface() const;
    inline void set_surface(Surface *surface);

    virtual Rgb first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler) const;
    virtual Rgb irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const;

protected:
    Surface *surface;
    Rgb radiance;
};

inline Surface *AreaLight::get_surface() const {
    return surface;
}

inline void AreaLight::set_surface(Surface *surface) {
    this->surface = surface;
}

inline bool AreaLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return surface->hit(r, tmin, tmax, time, record);
}

inline bool AreaLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return surface->hit(r, tmin, tmax, time);
}

inline Float AreaLight::area() const {
    return surface->area(); 
}

inline Rgb AreaLight::power() const {
    return radiance * area() * A_PI; // integrate L cos(theta) sin(theta) d(theta) d(phi) = pi
}

inline void AreaLight::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    surface->sample(rd, p, n, pdf);    
}

inline void AreaLight::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const {
    surface->sample(rd, p, n, pdf, patch);
}

inline void AreaLight::sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf) const {
    Float pdf_p;
    Float pdf_wo;
    this->sample(rd, p, n, wo, radiance, pdf_p, pdf_wo);
    pdf = pdf_p * pdf_wo;
}

inline Float AreaLight::pdf(const Vec3 &p, const Receiver &patch) const {
    return surface->pdf(p, patch);
}

inline Float AreaLight::pdf(const Vec3 &p) const {
    return surface->pdf(p);
}
    
inline Float AreaLight::pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo) const {
    if (dot(wo, n) < ZERO_EPSILON) return 0.0f;

    Float pdf_p = surface->pdf(p);
    Float pdf_wo;
    if (pdf_p <= 0.0f) 
        pdf_wo = 0.0f;
    else
        pdf_wo = dot(wo, n) * INV_PI;

    return pdf_p * pdf_wo;
}

inline void AreaLight::pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo, Float &pdf_p, Float &pdf_wo) const {
    if (dot(wo, n) < ZERO_EPSILON) {
        pdf_p = 0.0f;
        pdf_wo = 0.0f;
    } else {
        pdf_p = surface->pdf(p);
        if (pdf_p <= 0.0f)
            pdf_wo = 0.0f;
        else
            pdf_wo = dot(wo, n) * INV_PI;
    }
}

inline Float AreaLight::pdf_wo(const Vec3 &p, const Vec3 &wo) const {
    Vec3 n = surface->get_shape()->shading_normal(p);
    return std::max(dot(wo, n) * INV_PI, 0.0f);
}

} // end namespace Renzoku

#endif

