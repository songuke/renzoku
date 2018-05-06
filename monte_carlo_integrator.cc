#include "monte_carlo_integrator.h"
#include "scene.h"
#include "material.h"
#include "area_light.h"
#include "aggregate.h"
#include "image.h"
#include "frame.h"

namespace Renzoku {

MonteCarloIntegrator::MonteCarloIntegrator() {
    suffix = "mc";    
}

void MonteCarloIntegrator::initialize(Scene *scene) {    
    Integrator::initialize(scene);

    Size2 img_size = scene->get_image_size();
    target = new ImageFloat(0.0f, img_size.height, img_size.width, 3);
    weight = new ImageFloat(0.0f, img_size.height, img_size.width, 1);
}

Rgb MonteCarloIntegrator::radiance(const Ray &r) {
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    
    HitRecord rec;    
    if (! agg->hit(r, tmin, max_tmax, tick, rec)) return DefaultRgb::black;
    
    Vec3 p = r.org() + rec.t * r.dir();    
    Vec3 wo = unit_vector(-r.dir());
    
    Receiver recv(r, rec);
    return radiance(recv);
}

void MonteCarloIntegrator::radiance(PrimaryRay &pr) {    
    Rgb Lo = radiance(pr.r);
    pr.val = Lo;
}

void MonteCarloIntegrator::integrate(PrimaryRayBundle &b, FrameBuffer *frame) {    
    for (int i = 0; i < b.rays.size(); ++i) {
        radiance(b.rays[i]);        
    }
    
    Float filter_weight = 0.25f;    // 2x2 non-overlap box filter at each pixel
    // TODO: can retrieve subpixel and determine overlapping filter weight. 

    ImageFloat *img = frame->get_current_buffer();
    for (int k = 0; k < b.rays.size(); ++k) {
        int i = b.rays[k].s.y;
        int j = b.rays[k].s.x;
        Rgb v = b.rays[k].val;
        
        (*target)(i, j, 0) += filter_weight * v.red();
        (*target)(i, j, 1) += filter_weight * v.green();
        (*target)(i, j, 2) += filter_weight * v.blue();
        (*weight)(i, j, 0) += filter_weight;

        if ((*weight)(i, j, 0) == 0) continue; 
        Float inv_weight = 1.0f / (*weight)(i, j, 0);
        (*img)(i, j, 0) = (*target)(i, j, 0) * inv_weight;
        (*img)(i, j, 1) = (*target)(i, j, 1) * inv_weight;
        (*img)(i, j, 2) = (*target)(i, j, 2) * inv_weight;
    }
}

Rgb MonteCarloIntegrator::radiance(const Receiver &r) {
    return DefaultRgb::black;
}

} // end namespace