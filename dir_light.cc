#include "dir_light.h"
#include "scene.h"
#include "shading_geometry.h"

namespace Renzoku {

DirectionalLight::DirectionalLight(const Vec3 &wi, const Rgb &emission) : color(emission), wi(unit_vector(wi)) {
}

bool DirectionalLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return false;
}
   
bool DirectionalLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return false;
}
Float DirectionalLight::area() const {
    return approx_area;
}

Rgb DirectionalLight::power() const {
    return color * approx_area;
}
    
void DirectionalLight::query(Vec3 &wi, Rgb &radiance) const {
    wi = this->wi;
    radiance = this->color;
}

void DirectionalLight::set_approximate_area(Scene *scene) {    
    BoundingSphere bs = scene->get_bounding_sphere();
    this->approx_area = bs.rad * bs.rad * A_PI;   // the disk area given by the sphere equator
}


Rgb DirectionalLight::first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler) const {
    Aggregate *agg = scene->get_aggregate();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();

    if (r.light || agg->hit(Ray(r.p, wi), tmin, max_tmax, tick)) return DefaultRgb::black;

    ShadingGeometry sg(r);    
    Rgb radiance = color * sg.eval(wi) * sg.cosine(wi);
    return radiance;    
}

Rgb DirectionalLight::irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const {
    Aggregate *agg = scene->get_aggregate();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();

    if (agg->hit(Ray(dg.p, wi), tmin, max_tmax, tick)) {
        sr.is_valid = false;
        return DefaultRgb::black;
    }
  
    // use absolute value because we do not consider surface material for irradiance
    Rgb irradiance = color * fabs(dot(wi, dg.n));
    return irradiance;    
}

} // end namespace Renzoku