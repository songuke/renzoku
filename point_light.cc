#include "point_light.h"
#include "scene.h"
#include "shading_geometry.h"
#include "form_factor.h"

namespace Renzoku {

PointLight::PointLight(const Vec3 &_p, const Rgb _intensity) 
    : p(_p), intensity(_intensity) {}

Vec3 PointLight::org() const { return p; }

Rgb PointLight::radiance() const { return intensity; }

bool PointLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return false;
}

bool PointLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return false;
}

Float PointLight::area() const {
    return 0.;
}

Rgb PointLight::power() const {
    return 2 * A_PI * intensity;
}

void PointLight::sample(Random &rd, Vec3 &p, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const {
    // uniform sample hemisphere
    Vec2 uv;
    uv.random(rd);
    Float cos_theta = uv[0];
    Float sin_theta = sqrt(1.0f - cos_theta * cos_theta);
    Float phi = 2 * A_PI * uv[1];
    
    Float x = sin_theta * cos(phi);
    Float y = sin_theta * sin(phi);
    Float z = cos_theta; 

    wo = Vec3(x, y, z);
    p = this->p;
    radiance = this->intensity;
    pdf_p = 1.0f;
    pdf_wo = 1.0f / TWO_PI;
}

void PointLight::sample(Random &rd, Vec3 &p, Float &pdf) const {
    p = this->p;
    pdf = 1.0f;
}

Float PointLight::pdf(const Vec3 &p) const {
    return 0.0f;
}

Float PointLight::pdf(const Vec3 &p, const Vec3 &wo) const {
    return 0.0f;
}

void PointLight::query(const Vec3 &wo, Rgb &radiance) const {
    radiance = this->intensity;
}

Rgb PointLight::first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();    

    Vec3 q = this->org();
    if (r.light || agg->hit(r.p, q, tmin, tick)) return DefaultRgb::black;

    ShadingGeometry sg(r);
    Vec3 unit_qp = unit_vector(r.p - q);    
    Float G = FormFactor::form_factor_abs(q, unit_qp, r.p, r.shading_n);
    Rgb radiance = intensity * G * sg.eval(-unit_qp);
    return radiance;
}

Rgb PointLight::irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();    

    Vec3 q = this->org();
    if (agg->hit(dg.p, q, tmin, tick)) {
        sr.is_valid = false;
        return DefaultRgb::black;
    }
    
    Vec3 unit_qp = unit_vector(dg.p - q);
    sr.wo = unit_qp;

    // We use se absolute value for form factor because source is a light, and considered two-sided by default. 
    // Irradiance do not consider material at the receiver, hence the absolute at the receiver. So overall, use absolute.
    Float G = FormFactor::form_factor_abs(q, unit_qp, dg.p, dg.n);    
    Rgb irradiance = intensity * G;
    return irradiance;
}
 
} // end namespace Renzoku
