#ifndef _FIRST_BOUNCE_H_
#define _FIRST_BOUNCE_H_

#include "common.h"
#include "monte_carlo_integrator.h"

#include "random.h"
#include "material.h"
#include "light.h"
#include "area_light.h"
#include "aggregate.h"
#include "scene.h"

#include "local_geometry.h"
#include "shading_geometry.h"

#include "log.h"

namespace Renzoku {

/**
 * \brief FirstBounce integrator calculates direct illumination from light sources.
 *
 * BRDF evaluation takes care of ray under an opaque surface (return black and zero pdf). 
 * Therefore, form factor and cosine angle can always be regarded as positive. 
 *
 */ 
class FirstBounce : public MonteCarloIntegrator {
public:
    FirstBounce();

    virtual void initialize(Scene *scene);
    inline virtual Rgb radiance(const Receiver &);

    /**
     * Direction sampler can be passed in for MIS use.
     */
    inline static Rgb radiance(Scene *scene, const Receiver &r, 
                               DirectionSampler *dir_sampler = NULL); 

protected:
    //int num_light_source_samples;       /// for multiple importance sampling that combines light source and BRDF sampling
    //int num_brdf_samples;               /// for multiple importance sampling that combines light source and BRDF sampling
    //int num_pixel_samples;

    //Sampler *sampler;
};


inline Rgb FirstBounce::radiance(const Receiver &r) {    
    return radiance(scene, r);
}

inline Rgb FirstBounce::radiance(Scene *scene, const Receiver &r, 
                                 DirectionSampler *dir_sampler) {
    if (r.light) {
        Rgb Le;
        r.light->query(r.p, r.shading_n, r.wo, Le);
        return Le;
    }
    
    Lights &lights = *scene->get_lights();
    Random &rd     = *scene->get_random();
    
    Float light_pdf;        
    int k = scene->sample_light(rd(), light_pdf);
    if (k < 0) return DefaultRgb::black;        
    
    Rgb Lo = lights[k]->first_bounce(scene, r, dir_sampler) / light_pdf;
    return Lo;
}

} // end namespace Renzoku

#endif

