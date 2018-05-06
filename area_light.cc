#include "area_light.h"
#include "light.h"
#include "surface.h"
#include "scene.h"
#include "form_factor.h"
#include "mis.h"
#include "shading_geometry.h"
#include "hitrecord.h"

#include "log.h"
#include "stats.h"

namespace Renzoku {

AreaLight::AreaLight(Surface *surface, const Rgb &radiance)
    : surface(surface), radiance(radiance) {
}

void AreaLight::query(const Vec3 &p, const Vec3 &n, const Vec3 &wo, Rgb &radiance) const {
    Float pdf = this->pdf(p, n, wo);
    if (pdf <= 0.0f) 
        radiance = DefaultRgb::black;
    else
        radiance = this->radiance;
}

void AreaLight::sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const {    
    surface->sample(rd, p, n, pdf_p);
    
    // cosine weighted sample a hemisphere
    Vec2 sampler;
    sampler.random(rd);
    Float phi = TWO_PI * sampler.x();
    Float sin_theta = sqrt(sampler.y());
    Float x = sin_theta * cos(phi);
    Float y = sin_theta * sin(phi);
    Float z = sqrt(1.f - sin_theta * sin_theta); 

    // transform into world space
    Onb uvn = surface->get_basis(p);
    wo = unit_vector(uvn.local_to_world(Vec3(x, y, z)));
        
    radiance = this->radiance;

    pdf_wo = z * INV_PI;
} 

Rgb AreaLight::first_bounce(Scene *scene, const Receiver &r, 
                            DirectionSampler *dir_sampler) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Random &rd     = *scene->get_random();

    if (r.light) {        
        if (r.light == this) {

            Rgb Le;
            r.light->query(r.p, r.shading_n, r.wo, Le);
            return Le;

        } else {
            return DefaultRgb::black;       // no light between two luminaires
        }
    }

    Rgb Lo;
    Vec3 p = r.p;
    Vec3 pn = r.shading_n;
    ShadingGeometry sg(r);

    Material *m = r.m;
    bool is_singular;
    is_singular = (dir_sampler && dir_sampler->is_valid() && dir_sampler->is_singular()) || 
                  (m->get_bsdf()->is_singular());
    
    if (! is_singular) {

        // sample a point on the area light
        Vec3 q, qn;     
        Float pdf;

        // TODO: feed a Sampler to the sample function, or pass a 2D generated from the Sampler
        // but it needs to be another sampler            
        this->sample(rd, q, qn, pdf, r);

        if (pdf > 0.0f && ! agg->hit(p, q, tmin, tick)) {
            Rgb Le;
            Vec3 wi = unit_vector(q - p);            
            this->query(q, qn, -wi, Le);
                        
            Float G = FormFactor::form_factor_abs(p, pn, q, qn);

            if (G > 0.0f) {
                Rgb brdf = sg.eval(wi);

                // pdf as if the light sample is generated from the sampler
                Float dir_sampler_pdf;
                if (dir_sampler)
                    dir_sampler_pdf = dir_sampler->pdf(r, wi);
                else
                    dir_sampler_pdf = m->pdf(LocalGeometry(r), r.wo, wi);
                
                // convert the pdf into solid angle unit for MIS weight
                Float weight1 = Mis::balance(pdf * (q - p).squared_length() / fabs(dot(unit_vector(p - q), qn)), 
                                             dir_sampler_pdf);
                Lo += (weight1 * G / pdf) * Le * brdf;
            }
        }
    }
    
    // sample a direction using direction sampler and try to hit the light            
    Vec3 wi;
    Float dir_sampler_pdf;    
    if (dir_sampler)         
        dir_sampler->sample(rd, r, wi, dir_sampler_pdf);
    else
        m->sample(rd, LocalGeometry(r), r.wo, wi, dir_sampler_pdf);

    HitRecord rectwo;
    if (dir_sampler_pdf > 0.0f &&
        agg->hit(Ray(p, wi), tmin, max_tmax, tick, rectwo) &&         
        rectwo.light == this) {
                                    
        Vec3 q = p + rectwo.t * wi;
        Vec3 qn = rectwo.shading_normal;
        Rgb Le;
        this->query(q, qn, -wi, Le);
        
        Float weight2 = 1.0f;
        if (! is_singular) {            
            // convert to solid angle measure
            Float light_pdf = this->pdf(q, r) * (q - p).squared_length() / fabs(dot(unit_vector(p - q), qn));        
            weight2 = Mis::balance(dir_sampler_pdf, light_pdf);
        }
                        
        Rgb brdf = sg.eval(wi);
        Lo += (weight2 * sg.cosine(wi) / dir_sampler_pdf) * Le * brdf;
    }
    return Lo;
}

Rgb AreaLight::irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const {
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();
    Random &rd     = *scene->get_random();
    
    Vec3 p = dg.p;
    Vec3 pn = dg.n;
        
    Vec3 q, qn;     
    Float pdf;        
    this->sample(rd, q, qn, pdf);

    if (pdf <= 0.0f || agg->hit(p, q, tmin, tick)) {
        sr.is_valid = false;
        return DefaultRgb::black;
    }
    
    sr.wo = unit_vector(p - q);
    sr.p = q;
    sr.n = qn;
    sr.has_location = true;
    sr.has_normal = true;

    Rgb Le;
    this->query(q, qn, sr.wo, Le);                        
    Float G = FormFactor::form_factor_abs(p, pn, q, qn);
    Rgb irradiance = Le * G / pdf;
    return irradiance;
}

} // end namespace Renzoku
