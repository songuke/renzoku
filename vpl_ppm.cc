#include "vpl_ppm.h"

#include "rgb.h"
#include "vec2.h"
#include "vec3.h"

#include "random.h"
#include "camera.h"
#include "material.h"
#include "light.h"
#include "shape.h"
#include "surface.h"
#include "aggregate.h"

#include "scene.h"
#include "image.h"

#include "area_light.h"
#include "brdf_point_light.h"

#include "first_bounce.h"

#include "stats.h"
#include "log.h"

namespace Renzoku {

/**
 * When the number of photons is low, the caustics are very blurry and cannot be seen.
 * 
 * Increase the number of photons for sharper caustics.
 */
VplPpm::VplPpm() 
    : num_nearest_photons(64),
      use_sppm(true),
      num_emitted_photons(0),
      alpha(0.67f)
{
    
}

VplPpm::~VplPpm() {
}

void VplPpm::initialize(Scene* scene) {
    MonteCarloIntegrator::initialize(scene);

    if (vplgen == NULL)
        vplgen = new VplGenerator();

    // do not store direct VPLs because there is no wi.
    vplgen->store_direct_vpls(false);

    Stats stats;
    
    Log::info() << "Generating VPLs ... " << endn;
    stats.tic();
    vplgen->spawn(scene, num_start_particles,
                  this->get_max_bounce() + 1, all_vpls);
    stats.toc();
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    Log::info() << "Building kd-tree ... " << endn;
    stats.tic();
    kdtree = new KdTree(all_vpls, *scene->get_random());
    //kdtree = new KdTree(all_vpls);
    stats.toc();
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    this->observable.attach_observer(this);

    domain_stats = NULL;
    if (use_sppm) {
        Size2 size = scene->get_image_size();
        int num_pixels = size.height * size.width;
        domain_stats = new DomainStats[num_pixels];
    }
    //num_emitted_photons = all_vpls.size();
    num_emitted_photons = num_start_particles;
    // Normalize across all starting particles, because photons on the same subpath accounts for bounces. 
    // more bounces means more light accumulated to the image. 
    // If we normalize across total photons, it means when number of bounces increases, the image will get darker
    // due to higher bounces convey little energy.
    // (which is not expected because bounce energy should accumulate, not average out).

    // Also, when more bounces are traced, more higher bounce and less energy are emitted. Given the same disk,
    // the accumulated energy is going to be less and the image looks darker. 

    // To overcome the above problem, use Russian roulette to trace photons. 
    // To avoid increasing too much variance, only use RR after 3 or 4 bounces. By assuming that higher bounce convey less energy
    // variance when increased will be less visible.

}

Rgb VplPpm::photon_gather(const Receiver &r, int eye_bounce, int max_bounce) {
    // receiver is not on light because we already check for that

    //
    // From Jensen's photon mapping course at SIGGRAPH 2008:
    // caustics         : gather directly at the hitpoint
    // interreflection  : sample a ray to find the next hitpoint and gather there (final gathering)
    // 
    // Here we simply render everything at the first hit point.
    //

    // for each hitpoint, find its k nearest neibour photons in the photon map    
    LightParticleHeap heap(r.p);
    kdtree->find_nearest(r.p, num_nearest_photons, heap);
    Float radius = (heap.top()->org() - r.p).length();
    
    Rgb Lo;
    Float inv_dA = INV_PI / (radius * radius);
    
    // cone filter to weigh photons near to p more for sharper edges
    Float k = 1.f;
    Float total_weight = 0.f;
    
    LocalGeometry dg(r);

    LightParticlePtrs *nearest = reinterpret_cast<LightParticlePtrs *>(&heap);
    for (int i = 0; i < nearest->size(); i++) {
        BrdfPointLight *vpl = (*nearest)[i]->get_brdf_point_light();
        
        if (eye_bounce + vpl->get_bounce() - 1 > max_bounce) continue;  // VPL at area light (0 bounce) is not stored because no wi.
        if (eye_bounce + vpl->get_bounce() - 1 == 1 && scene->get_integrator()->is_direct_lighting() == false) continue;
        if (vpl->get_bounce() == 0 && eye_bounce > 1) continue;         // avoid double count

        Float weight = std::max(1e-2f, 1.f - ((*nearest)[i]->org() - r.p).length() / (k * radius));                
        total_weight += weight;

        Rgb flux = vpl->power();
        Rgb brdf = r.m->eval(dg, r.wo, vpl->get_wi());
        Lo += brdf * flux * weight;
    }
    Lo *= inv_dA / total_weight / all_vpls.size();
        
    return Lo;
}

void VplPpm::on_complete_frame(int tid, ImageFloat *img) {
    if (tid != 0) return;

    Stats stats;

    Log::info() << "Cleaning up VPLs ... " << endn;
    stats.tic();
    all_vpls.clear();
    stats.toc();    
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    Log::info() << "Cleaning up kd-tree ... " << endn;
    stats.tic();
    delete kdtree; 
    stats.toc();
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    Log::info() << "Generating VPLs ... " << endn;
    stats.tic();
    // for photons, due to possible vertex merge with the first primary hit point in gathering, 
    // we need to shoot one more bounce.
    vplgen->spawn(scene, num_start_particles,
                  this->get_max_bounce() + 1, all_vpls);
    stats.toc();
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    Log::info() << "Building kd-tree ... " << endn;
    stats.tic();
    kdtree = new KdTree(all_vpls, *scene->get_random());
    //kdtree = new KdTree(all_vpls);
    stats.toc();
    Log::info() << "\t" << stats.elapsed() << " ms." << endn;

    //num_emitted_photons += all_vpls.size();
    num_emitted_photons += num_start_particles;
}

void VplPpm::radiance(PrimaryRay &pr) {
    Rgb Lo = mixed_gather(pr);
    pr.val = Lo;
}

Rgb VplPpm::mixed_gather(const PrimaryRay &pr) {
    
    Aggregate *agg = scene->get_aggregate();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Random &rd     = *scene->get_random();

    Rgb Lo = DefaultRgb::black;
    Rgb T = DefaultRgb::white;
    int bounce = 0;
    Ray r = pr.r;

    do {
        bounce++;
    
        HitRecord rec;
        if (! agg->hit(r, tmin, max_tmax, tick, rec)) break;
    
        LocalGeometry dg(r, rec);
        Vec3 p = dg.p;
        Vec3 pn = dg.n;
        Vec3 wo = unit_vector(-r.dir());

        if (rec.light) {
            if (bounce == 1 && scene->get_integrator()->is_direct_lighting()) {
                // since FirstBounce accounts for direct illumination (path length 2), 
                // here we should only query of path is of length 1 
                Rgb Le;
                rec.light->query(p, pn, -r.dir(), Le);
                Lo += T * Le;
                break;
            } 
        }

        Material *m = rec.material;
        if (m == NULL) break;
        Receiver recv(r, rec); 

        // photons merge vertex, so it already includes direct illumination

        // only gather directly if the surface is non-singular
        //if (m->get_bsdf()->is_singular() == false) {        
        if (m->get_bsdf()->get_bsdf_type() == Bsdf::LAMBERTIAN) {
            if (use_sppm) {
                int width = scene->get_image_size().width;
                int pixel_index = pr.s.y * width + pr.s.x;
                Lo += T * photon_gather(recv, pixel_index, bounce, this->get_max_bounce());
            } else {
                Lo += T * photon_gather(recv, bounce, this->get_max_bounce());
            }
            break;
        }
    
        // otherwise, continue to tracing until a non-specular vertex is found
        Lo += T * FirstBounce::radiance(scene, recv);

        Vec3 wi;
        Float pdf_wi;
        Rgb brdf = m->sample(rd, dg, -r.dir(), wi, pdf_wi);
        if (pdf_wi == 0.0f || brdf == DefaultRgb::black) break;
    
        T *= fabs(dot(pn, wi)) * brdf / pdf_wi;
        r = Ray(p, wi);

    } while (bounce < max_bounce);
    
    return Lo;
}

Rgb VplPpm::photon_gather(const Receiver &r, int pixel_index, int eye_bounce, int max_bounce) {
    // receiver is not on light as we already check for this before
    // if (r.m == NULL) return DefaultRgb::black;
        
    LightParticleHeap heap(r.p);
    DomainStats &ds = domain_stats[pixel_index];
    
    if (ds.radius < 0.0f) {  // first round            
        kdtree->find_nearest(r.p, num_nearest_photons, heap);
        if (heap.size() > 0)
            ds.radius = (heap.top()->org() - r.p).length();
        else
            ds.radius = scene->get_bounding_box().diagonal() * 0.01f;
    }
    kdtree->find_nearest(r.p, ds.radius, heap);
    
    LocalGeometry dg(r);

    Rgb new_radiance;
    LightParticlePtrs *nearest = reinterpret_cast<LightParticlePtrs *>(&heap);
    int new_photons = 0;
    for (int i = 0; i < nearest->size(); i++) {
        BrdfPointLight *vpl = (*nearest)[i]->get_brdf_point_light();
        
        if (eye_bounce + vpl->get_bounce() - 1 > max_bounce) continue;  // VPL at area light (0 bounce) is not stored because no wi.
        if (eye_bounce + vpl->get_bounce() - 1 == 1 && scene->get_integrator()->is_direct_lighting() == false) continue;        
        if (vpl->get_bounce() == 0 && eye_bounce > 1) continue;         // avoid double count

        Rgb flux = vpl->power();
        Rgb brdf = r.m->eval(dg, r.wo, vpl->get_wi());
        new_radiance += brdf * flux;

        new_photons++;
    }
        
    /**
     * Update domain stats
     */    
    if (new_photons > 0) {
        Float ratio = (ds.num_photons + alpha * new_photons) / (ds.num_photons + new_photons);
        ds.radius *= sqrt(ratio);
        ds.num_photons += alpha * new_photons;
        ds.acc_radiance = (ds.acc_radiance + new_radiance) * ratio;
    }
    Float inv_dA = INV_PI / (ds.radius * ds.radius);
    Rgb Lo = ds.acc_radiance * inv_dA / num_emitted_photons;

    return Lo;
}



} // end namespace
