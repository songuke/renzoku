#include "path_tracing.h"

#include "rgb.h"
#include "vec2.h"
#include "vec3.h"

#include "random.h"
#include "camera.h"
#include "material.h"
#include "light.h"
#include "area_light.h"
#include "dir_light.h"
#include "point_light.h"
#include "spot_light.h"
#include "shape.h"
#include "surface.h"
#include "aggregate.h"

#include "scene.h"
#include "image.h"
#include "sampler.h"

#include "monte_carlo_integrator.h"

#include "form_factor.h"
#include "first_bounce.h"
#include "shading_geometry.h"

#include "frame.h"

#include "log.h"
#include "stats.h"

namespace Renzoku {

PathTracing::PathTracing() { 
    suffix = "pt";
}

void PathTracing::initialize(Scene *scene) {
    MonteCarloIntegrator::initialize(scene);
}

/**
 * Path tracing with returning radiance value at each bounce.
 * 
 * Written at 4 Allenby, Wellington. 
 */
static Rgb trace_wellington(Scene *scene,
                            const Receiver &r, int bounce, int max_bounce) {
    
    Lights &lights = *scene->get_lights();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Integrator* integrator = scene->get_integrator();

    Vec3 p = r.p;
    Vec3 pn = r.shading_n;
    Vec3 wo = r.wo;

    // special handling for maximum zero bounce
    if (max_bounce == 0) {
        if (r.light && integrator->is_direct_lighting()) { 
            Rgb Le;
            r.light->query(p, pn, wo, Le);
            return Le;
        }
        return DefaultRgb::black;        
    }

    // other cases
    if (bounce > max_bounce) return DefaultRgb::black;

    if (r.light) { // ray hits the light directly        
        if (bounce == 1 && integrator->is_direct_lighting()) {
            Rgb Le;
            r.light->query(p, pn, wo, Le);
            return Le;
        } else {
            // ignore double count of direct illumination
            return DefaultRgb::black;
        }
    }

    // ray hits a surface
    Material *m = r.m;
    LocalGeometry dg(r);
    ShadingGeometry sg(r);
    
    Rgb radiance;
    if (bounce > 1 || (bounce == 1 && integrator->is_direct_lighting())) 
        radiance = FirstBounce::radiance(scene, r);
        
    if (bounce >= max_bounce) return radiance;        // early terminate if maximum path length is reached

    HitRecord rec;    
    Vec3 wi;
    Float pdf_wi;
    Rgb brdf = sg.sample(rd, wi, pdf_wi);
    if (pdf_wi == 0.0f || brdf == DefaultRgb::black) return radiance;

    Ray ray(p, wi);
    if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
        if (rec.light)  // avoid double count
            return radiance;

        Receiver recv(ray, rec);
        Rgb Li = trace_wellington(scene, recv, bounce + 1, max_bounce);

        radiance += Li * brdf * sg.cosine(wi) / pdf_wi;
    }
    return radiance;
}

Rgb PathTracing::trace(Scene *scene, const Receiver &r, int max_bounce) {
    return trace_wellington(scene, r, 1, max_bounce);
}

Rgb PathTracing::radiance(const Receiver &r) {
    Rgb Lo = trace_wellington(scene, r, 1, this->get_max_bounce());    
    return Lo;
}

static int frame_index = 0;
void PathTracing::on_frame_end() {    
    frame_index++;
    if (frame_index == 1 || (frame_index % 5 == 0)) {
        ostringstream oss;
        oss << scene->get_output_folder() << "/" << scene->get_name() << "_sample" << frame_index << ".exr";
        ImageFloat img(*scene->get_frame_buffer()->get_last_complete_frame_buffer());
    
        string file = oss.str();
        if (img.save(file.c_str()))
            Log::info() << file << " saved." << endn;
    }
}

} // end namespace Renzoku
