#include "vpl_generator.h"

#include "scene.h"
#include "camera.h"

#include "dir_light.h"
#include "area_light.h"
#include "spot_light.h"
#include "point_light.h"

#include "random.h"
#include "aggregate.h"
#include "hitrecord.h"
#include "integrator.h"
#include "path.h"

#include "shading_geometry.h"

#include "stats.h"
#include "log.h"

namespace Renzoku {

VplGenerator::VplGenerator() {
}

/**
 * Handle singularity in light source:
 * 
 * 1. Directional light has a delta distribution of direction. 
 * In BDPT, a path vertex can be created to represent the directional light to take care of path weight,
 * but connection to this path vertex is not valid. 
 * Direct illumination due to directional light can be evaluated using other techniques such as unidirectional (PT, LT).
 * Since connection is invalid, in VPL, we don't create this path vertex.
 *
 * 2. Environment light.
 * We can represent points on a large sphere as VPLs on the environment light. Such VPLs is connectable.
 *
 * For consistency and a unified shader, we only use VPL to evaluate indirect illumination, where surface point and direction are well defined. 
 * All direct illumination can be evaluated separately.
 */

void VplGenerator::trace_vpl(Scene *scene, const Ray &r, Rgb T, int bounce, int max_bounce, 
                             BrdfPointLights &vpls, const BrdfPointLightRef &prev, Rgb path_throughput, Float path_pdf) {
    // NOTE: trace 1 bounce less than path tracing
    if (bounce >= max_bounce) return;

    Lights &lights = *scene->get_lights();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();

    HitRecord rec;    
    if (! agg->hit(r, tmin, max_tmax, tick, rec)) return;
    
    if (rec.light) return;

    LocalGeometry dg(r, rec);
    Vec3 q = r.org();
    Vec3 p = dg.p;
    Vec3 pn = dg.n;
    Material *m = rec.material;    
    Vec3 wi = unit_vector(-r.dir());

    BrdfPointLight pl(dg, T, m, wi, prev);
    pl.set_bounce(bounce);
    pl.set_shape(rec.shape);
	pl.set_distance(rec.t);
    pl.set_throughput(path_throughput);
    pl.set_pdf(path_pdf);

    /*
    // LSD* subpath only
    if (bounce == 1 && m->get_bsdf_type() == Material::LAMBERTIAN) return;
    if (m->get_bsdf_type() == Material::LAMBERTIAN)
        vpls.push_back(pl);
    */
    vpls.push_back(pl);
    
    // for bounces higher than 4, start russian roulette
    Float continuation = 1.0f;
    if (bounce >= 4) {
        continuation = std::min(m->get_representative_color().avg(), 0.99f);

        if (rd() > continuation) {
            return;            
        } 
    }

    Vec3 wo;
    Float pdf;
    Rgb brdf = m->sample(rd, dg, wi, wo, pdf);
    
	if (pdf <= 0.0f || brdf == DefaultRgb::black) return;

    pdf *= continuation;
    
    // consistent geometry normal check
    //if (dot(rec.normal, wo) <= ZERO_EPSILON) return;
                
    // prepare T with outgoing cosine at q    
	T *= brdf * fabs(dot(pn, wo)) / pdf; // very much like tracing potential (in light tracing)
    
    path_throughput *= brdf * fabs(dot(pn, wo));
    path_pdf *= pdf;
	trace_vpl(scene, Ray(p, wo), T, bounce + 1, max_bounce, vpls, BrdfPointLightRef(vpls.size() - 1), path_throughput, path_pdf);
    return;
    
	// TODO: use adaptive diffuse component for each surface
    /*
    Float absorption = rd();
	Float average_albedo = std::max(0.5f, scene->get_average_albedo());
	if (absorption < average_albedo) { // continue
		pdf *= 1.0f / (absorption);

		// prepare T with outgoing cosine at q    
		T *= brdf * dot(pn, wo) / pdf; // very much like tracing potential (in light tracing)
		trace_vpl(scene, Ray(p, wo), T, bounce + 1, max_bounce, vpls, pl, path_throughput, path_pdf);

	} else {
		return;
    }*/
}
    
void VplGenerator::generate_vpl(Scene *scene, int num_light_samples,
                                int max_bounce, BrdfPointLights &vpls, bool is_photon) {

    Random &rd = *scene->get_random();

    for (int i = 0; i < num_light_samples; ++i) {
        PathNode vertex;

        // TODO: extra parameters such as shape, path_throughput, path_pdf might be removed in the future.
        DirectVertexSamplingRecord sr;
        if (! SubpathGenerator::sample_direct_vertex(scene, vertex, sr)) continue;

        // create bounce-0 VPL to represent direct lighting
        if (sr.emitter_dg_valid) {
            BrdfPointLight emitter_vpl(sr.emitter_dg, sr.emitter_power, sr.emitter);
            emitter_vpl.set_bounce(0);
            //vpls.push_back(emitter_vpl);    // Apr 06: since we don't need to guide direct lighting, no bounce-0 VPL.
        }

        if (vertex.pdf_dA == 0.0f) continue;    // invalid vertex
        
        // turn the vertex into a VPL, and continue tracing        
        Rgb throughput = sr.vertex1_throughput;
        Float path_pdf = sr.vertex1_pdf_dA;

        if (max_bounce <= 1) continue;

        // for photon, this will approximate direct illumination
        // we want to skip this if only indirect illumination is considered
        if (! is_photon) {
            BrdfPointLight pl(vertex.dg, vertex.light_contrib(), vertex.material, vertex.actual_wi, sr.emitter);        
            pl.set_bounce(1);
            pl.set_shape(sr.vertex1_shape);
		    pl.set_distance((vertex.dg.p - sr.emitter_dg.p).length());
            pl.set_throughput(throughput);
            pl.set_pdf(path_pdf);        
            vpls.push_back(pl);
        }

        Material *m = vertex.material;
        Vec3 wi;
        Float pdf_wi;
        Rgb brdf = m->sample(rd, vertex.dg, vertex.actual_wi, wi, pdf_wi);
        if (pdf_wi <= 0.0f || brdf == DefaultRgb::black) continue;

        ShadingGeometry sg(vertex.actual_wi, vertex.material, vertex.dg);

        throughput *= brdf * sg.cosine(wi);
        path_pdf *= pdf_wi;
        Rgb contrib = vertex.light_contrib() * brdf * sg.cosine(wi) / pdf_wi;

        trace_vpl(scene, Ray(vertex.dg.p, wi), contrib, 2, max_bounce, vpls, BrdfPointLightRef(vpls.size() - 1), throughput, path_pdf);
    }

}

void VplGenerator::spawn(Scene *scene, int num_start_particles, 
                         int max_bounce,
                         BrdfPointLights &all_vpls, bool is_photon) {
    Stats stats;
    stats.tic();    

    int num_light_samples = num_start_particles;

    Log::info() << "Generating VPLs: " << num_light_samples << " samples..." << endn;    
    generate_vpl(scene, num_light_samples, max_bounce, all_vpls, is_photon);
    for (int i = 0; i < all_vpls.size(); ++i) {
        if (is_store_direct_vpls) {
            if (all_vpls[i].get_bounce() == 0)
                all_vpls[i].set_tag(VPL_TAG_NEW_GROUP);
        } else {            
            if (all_vpls[i].get_bounce() == 1)
                all_vpls[i].set_tag(VPL_TAG_NEW_GROUP);        
        }
    }

    Log::info() << "Total VPLs                        : " << all_vpls.size() << endn;
}

} // end namespace
