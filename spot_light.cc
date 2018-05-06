#include "spot_light.h"
#include "scene.h"
#include "form_factor.h"
#include "shading_geometry.h"
#include "log.h"

namespace Renzoku {

SpotLight::SpotLight(const Vec3 &_p, const Vec3 &_n, Float phi, Float theta, Float falloff, const Rgb &_intensity) 
    : pos(_p), normal(_n), phi(phi), theta(theta), falloff(falloff), intensity(_intensity) 
{
    onb.init_from_n(normal);
        
    cos_half_theta = cos(theta / 2);
    cos_half_phi = cos(phi / 2);    
    cone_size = TWO_PI * (1.0f - cos_half_phi);
    diff_cos = cos_half_theta - cos_half_phi;
}

Vec3 SpotLight::org() const { return pos; }
Vec3 SpotLight::direction() const { return normal; }

bool SpotLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return false;
}

bool SpotLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return false;
}

Float SpotLight::area() const {
    return 0.;
}

Rgb SpotLight::power() const {
    // TODO: this is just the upper bound. 
    // need to integrate with fall off
    return cone_size * intensity;
}

void SpotLight::sample(Random &rd, Vec3 &p, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const {
    // uniform sample the cone by phi angle
    Vec2 delta;
    delta.random(rd);

    Float cos_theta = 1.0f - delta.x() * (1.0f - cos_half_phi);
    Float sin_theta = sqrt(1.0f - cos_theta * cos_theta);
    Float phi = TWO_PI * delta.y();
    
    wo = onb.local_to_world(Vec3(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta));    
    p = this->pos;    
    pdf_p = 1.0f;
    pdf_wo = 1.0f / cone_size;
    this->query(wo, radiance);
}

void SpotLight::sample(Random &rd, Vec3 &p, Float &pdf) const {
    p = this->pos;
    pdf = 1.0f;
}

Float SpotLight::pdf(const Vec3 &p) const {
    return 0.0f;
}

Float SpotLight::pdf(const Vec3 &p, const Vec3 &wo) const {
    return 0.0f;
}

void SpotLight::query(const Vec3 &wo, Rgb &radiance) const {    
    if (dot(wo, normal) >= cos_half_theta) {
        radiance = this->intensity;
    } else if (dot(wo, normal) > cos_half_phi) {            
        Float ratio = (dot(wo, normal) - cos_half_phi) / diff_cos;
        Float factor = pow(ratio, falloff);
        radiance = this->intensity * factor;
    } else {
        radiance = DefaultRgb::black;
    }    
}

Rgb SpotLight::first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();  

    Vec3 q = this->org();
    if (r.light || agg->hit(r.p, q, tmin, tick)) return DefaultRgb::black;

    Vec3 unit_qp = unit_vector(r.p - q);
    Rgb Le;
	this->query(unit_qp, Le);
        
    ShadingGeometry sg(r);
    Float G = FormFactor::form_factor_abs(q, unit_qp, r.p, r.shading_n);
    Rgb radiance = Le * G * sg.eval(-unit_qp);
    return radiance;
}

Rgb SpotLight::irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();  

    Vec3 q = this->org();
    if (agg->hit(dg.p, q, tmin, tick)) {
        sr.is_valid = false;
        return DefaultRgb::black;
    }

    Vec3 unit_qp = unit_vector(dg.p - q);
    Rgb Le;
	this->query(unit_qp, Le);
    
    sr.wo = unit_qp;

    Float G = FormFactor::form_factor_abs(q, unit_qp, dg.p, dg.n);
    Rgb irradiance = Le * G;
    return irradiance;
}

} // end namespace Renzoku
