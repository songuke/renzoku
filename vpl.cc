#include "vpl.h"

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
#include "frame.h"

#include "area_light.h"
#include "brdf_point_light.h"
#include "dir_light.h"

#include "form_factor.h"
#include "first_bounce.h"

#include "stats.h"
#include "log.h"

namespace Renzoku {


static void visualize_dA(Scene *scene, const Ray &r, Rgb &radiance);

void VirtualPointLight::initialize(Scene *scene) {
    MonteCarloIntegrator::initialize(scene);

    Size2 size = scene->get_camera()->get_image_size();
    img_bias = new ImageFloat(0.f, size.height, size.width, 3);
    img_vpl  = new ImageFloat(0.f, size.height, size.width, 3);

    frame_index = 0;

    iteration = 0;
    num_evaluated_vpls = 0;
    num_evaluated_bias_samples = 1;    // for the first round
    
    num_emitted_photons = 0;    
    total_bounce1 = 0;

    spawn_vpls();

    switch (mis_type) {
        case METRO_ONLY:
            scene->set_output_folder("data/metro_only");
        break;
        case METRO_RANDOM:
            scene->set_output_folder("data/metro_random");
        break;
        case VORBA_ONLY:
            scene->set_output_folder("data/vorba_only");
        break;
        case VORBA_RANDOM:
            scene->set_output_folder("data/vorba_random");
        break;
    }
    Log::info() << "Scene output folder: " << scene->get_output_folder() << endn;
}

void VirtualPointLight::spawn_vpls() {
    Stats stats;

    stats.tic();
    if (vplgen == NULL)
        vplgen = new VplGenerator();

    vplgen->store_direct_vpls(false);

    if (estimate_mode == VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE) {
        // when use as photon, we need to trace 1 more bounce as compared to VPL
        // and skip storing the first bounce photon as it approximates direct illumination
        vplgen->spawn(scene, num_start_particles,
            this->get_max_bounce() + 1, all_vpls, true);

    } else {
        vplgen->spawn(scene, num_start_particles, 
                  this->get_max_bounce(), all_vpls);
    }

    if (all_vpls.size() <= 0) {
        Log::critical("No VPLs generated.", ExitCode::EXIT_CODE_SCENE_FILE_ERROR);
    }

    stats.toc();
    stats.print_elapsed_milliseconds();
        
    num_emitted_photons += num_start_particles;
    for (int i = 0; i < all_vpls.size(); ++i) {
        if (all_vpls[i].get_bounce() == 1)
            total_bounce1++;
    }

    switch (estimate_mode) {
        case VplRadianceEstimate::LIGHTCUTS:  
            Log::error() << "Lightcuts not supported." << endn;
            /*   
            Log::info() << "Initializing lightcuts ... " << endn;
            stats.tic();
            lightcuts.initialize(scene, all_vpls);
            lightcuts.set_clamp_threshold(this->clamp_threshold);
            stats.toc();
            Log::info() << "Initializing lightcuts done." << endn;
            stats.print_elapsed_seconds();
            */
            break;

        case VplRadianceEstimate::MRCS:
            Log::info() << "Initializing MRCS ... " << endn;
            stats.tic();        
            mrcs.initialize(scene, this, all_vpls);
            stats.toc();
            Log::info() << "Initializing MRCS done." << endn;
            stats.print_elapsed_seconds();
            break;

        case VplRadianceEstimate::LIGHTSLICE:  
        case VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE:          
            Log::info() << "Initializing LightSlice ... " << endn;
            stats.tic();        
            lightslice.initialize(scene, this, all_vpls);
            stats.toc();
            Log::info() << "Initializing LightSlice done." << endn;
            stats.print_elapsed_seconds();
            break;

        case VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE:
            Log::info() << "Initializing LibImportance ... " << endn;
            stats.tic();
			importance_manager = new Vorba::ImportanceManager(scene, *scene->get_param_dict());
            importance_manager->initialize(all_vpls);
            stats.toc();
            stats.print_elapsed_seconds();
            break;

        case VplRadianceEstimate::POWER_SAMPLING: {
            vector<Float> power_arr(all_vpls.size());
            for (int i = 0; i < all_vpls.size(); ++i) {
                power_arr[i] = all_vpls[i].power().luminance();
            }
            powers.set_distribution(power_arr);
            
            Random &rd = *scene->get_random();
            powers.sample(rd, light_index, pdf_light);        
                    
            Log::info() << "Begin light index: " << light_index << " " << pdf_light << endn;
            num_power_samples_so_far = 1;

            while (pdf_light <= 0.0f) {
                powers.sample(rd, light_index, pdf_light);     
                Log::error() << "Resampling..." << endn;
            }

            break;
        }

        case VplRadianceEstimate::NONE:
            light_index = 0;
            break;
        
        default:
            Log::warn() << "Failed to initialize estimate mode " << estimate_mode << endn;
            break;
    }

    standard_offset    = 0;                     // start gathering from VPL 0    
       
    next_batch         = false;

}

VirtualPointLight::VirtualPointLight() {
    suffix = "vpl";

    this->clamp_mode = VplRadianceClamp::FORM_FACTOR;
    this->is_bias_compensation  = false;

    // default values for clamping
    this->set_clamping_threshold(1e-3f);    
    this->set_incident_radiance_clamping_threshold(1e-3f);

    clamp_relax                     = 0.66f;    
    iteration                       = 0;
    num_bias_compensation_samples   = 256;

    img_vpl  = NULL;
    img_bias = NULL;

    estimate_mode = VplRadianceEstimate::POWER_SAMPLING;

    num_virtual_indirect_samples = 16;

	importance_manager = NULL;
}

VirtualPointLight::~VirtualPointLight() {
    if (img_bias) 
        delete img_bias;
    if (img_vpl)
        delete img_vpl;

	if (importance_manager)
		delete importance_manager;
}

static void clamp_term(RadianceTerms &term, Float threshold, VplRadianceClamp::Type mode) {
    switch (mode) {
        case VplRadianceClamp::REFLECTIVITY:
            term.eval_clamp_reflectivity(threshold);
            break;
        case VplRadianceClamp::RADIANCE:
            term.eval_clamp_radiance(threshold);
            break;
        case VplRadianceClamp::FORM_FACTOR:
            term.eval_clamp_form_factor(threshold);
            break;
        case VplRadianceClamp::NONE:            
        default:            
            term.eval_no_clamp();
            break;
    }
}

/**
 * Estimate the bias energy by path tracing from the surface gathering point.
 */
Rgb VirtualPointLight::compute_bias_energy(Scene *scene, 
                                           const LocalGeometry &dg, 
                                           const Vec3 &wo, Material *m) const {

    Lights &lights = *scene->get_lights();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();

    Rgb bias_energy;
    
    bool terminated = false;        

    Vec3 x       = dg.p;
    Vec3 xn      = dg.n;
    Vec3 wo_x    = wo;
    Material *mx = m;
    Onb uvn_x    = dg.uvn;

    Rgb throughput = DefaultRgb::white;
    int bounce = 1;
    while (! terminated) {        
        if (bounce >= this->get_max_bounce()) { // terminate 
            Vec2 sample; sample.random(rd);
            Float light_pdf;
            int k = scene->sample_light(sample.x(), light_pdf);
            if (k < 0) {
                throughput = DefaultRgb::black;
                terminated = true;
                break;
            }

            Light *light = lights[k];

            if (light->get_light_type() == Light::AREA_LIGHT) {
                AreaLight *area_light = (AreaLight *)light;
                Vec3 q, qn, wo_q;
                Float pdf_q;
                Rgb Le;                
                area_light->sample(rd, q, qn, pdf_q);
                
                if (agg->hit(x, q, tmin, tick)) { // blocked
                    throughput = Rgb();
                    terminated = true;
                    break;
                }
                
                wo_q = unit_vector(x - q);
                area_light->query(q, qn, wo_q, Le);

                Float G = FormFactor::form_factor_abs(x, xn, q, qn);             
                //Float bias_factor = std::max(G - clamp_factor, 0.f) / G;                

                Rgb brdf_x = mx->eval(dg, wo_x, -wo_q);

                //throughput *= Le * brdf_x * dot(xn, -wo_q) * bias_factor / pdf_q;
                throughput *= Le * brdf_x * G / (light_pdf * pdf_q);
                terminated = true;
                break;
            }
        }
            
        // otherwise, sample a direction based on BRDF and continue
        Vec3 wi;
        Float pdf;
        Rgb brdf = mx->sample(rd, dg, wo_x, wi, pdf);        

        Ray r(x, wi);
        HitRecord rec;
            
        /**
        * The max operator when computing the bias factor indicates that 
        * when the form factor G < clamp_factor, the bias compensation energy is zero
        * (or no bias has occured previously). 
        */

        if (! agg->hit(r, tmin, max_tmax, tick, rec)) {
            throughput = DefaultRgb::black;
            terminated = true;
            break;
        }

        if (rec.light) {    // hit light
            throughput = DefaultRgb::black;
            /*
            Vec3 p = r.org() + rec.t * r.dir();
            // get light radiance
            Rgb Le;
            rec.light->query(p, -r.dir(), Le);
                    
            Float G = form_factor(x, xn, p, rec.normal);
            if (G <= 0) {
                throughput = Rgb();
                terminated = true;
                break;
            }

            if (bounce == 1) {                
                Float bias_factor = std::max(G - clamp_threshold, 0.f) / G;

                throughput *= Le * brdf * dot(xn, wi) * bias_factor / pdf;
            } else {
                throughput *= Le * brdf * dot(xn, wi) / pdf;
            }*/
            terminated = true;
            break;
        }                       
        if (terminated) break;
                        
        Vec3 y = x + rec.t * wi;
        Vec3 yn = rec.normal;

        Float G = FormFactor::form_factor_abs(x, xn, y, yn);
        
        if (bounce == 1) { // do bias compensation from first hit from gathering surface only            
            Float bias_factor = std::max(G - clamp_threshold, 0.f) / G;
            throughput *= brdf * dot(xn, wi) * bias_factor / pdf;
        } else { // do normal path tracing
            throughput *= brdf * dot(xn, wi) / pdf;
        }

        // terminate early if throughput becomes zero
        if (throughput == DefaultRgb::black) break;
            
        // next path segment or terminate
        ++bounce;                         
        x = y;
        xn = yn;
        mx = rec.material;
        wo_x = -wi;
        uvn_x = rec.uvn;            
    }        
    bias_energy += throughput;
    
    return bias_energy;
}

static int sample_vpl(Random &rd, Float *cdf, int n) {
    Float u = rd();
    for (int i = 0; i < n; ++i) {
        if (u >= cdf[i] && u < cdf[i + 1]) {
            return i;
        }
    }
    cout << "Last cdf: " << cdf[n] << endl;
    cout << "U: " << u << endl;
    cout << "Found 1 in random. Exiting..." << endl;
    exit(1);
    return -1;
}

void VirtualPointLight::debug_vpl(BrdfPointLights &vpls, const char *file) {
    // write all VPL points to a PLY file.
    FILE *f = fopen(file, "w");
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "element vertex %d\n", vpls.size());
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "end_header\n");
    for (int i = 0; i < vpls.size(); ++i) {
        Vec3 p = vpls[i].org();
        fprintf(f, "%f %f %f\n", p.x(), p.y(), p.z());
    }
    fclose(f); 
}

Rgb VirtualPointLight::radiance(Scene *scene, const Receiver &r, const BrdfPointLight &light, bool &visible) {
    LocalGeometry dg(r);
    Rgb Lo = gather_point(scene, dg, r.wo, r.m, &light, visible);
    return Lo;
}

Rgb VirtualPointLight::incident_radiance(Scene *scene, const Receiver &r, const BrdfPointLight &vpl) {
    Vec3 p = r.p;
    Vec3 pn = r.n;
    Vec3 q = vpl.org();
    Vec3 qn = vpl.normal();        
    
    Float G = FormFactor::form_factor_abs(p, pn, q, qn);    
    // no need to leave out the cosine because we don't use RIS
    //
    //Vec3 pq = q - p;    
    //Float G = fabs(dot(qn, unit_vector(pq))) / dot(pq, pq);


    Vec3 wi = unit_vector(q - p);
    
    Rgb light_brdf = vpl.brdf(-wi);    
    Rgb light_power = vpl.power();

    // no clamp
    //return light_power * light_brdf * G;

    RadianceTerms term(light_power, DefaultRgb::white, light_brdf, G);
    
    // since we use the incident radiance for sampling, it is importance to clamp away some singular peaks
    // otherwise it will cause noise
        
    // the threshold for reflectivity and radiance clamp needs to be adjusted to account for the lack of the gather point BRDF
    clamp_term(term, incident_radiance_clamp_threshold, clamp_mode);
    return term.radiance;
}

void VirtualPointLight::radiance(PrimaryRay &pr) {
    Rgb Lo;

    switch (estimate_mode) {
        case VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE:
        case VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE:
            Lo = guided_path_trace(pr);
            break;

        case VplRadianceEstimate::BIAS_ONLY:
            Lo = gather_bias(pr);
            break;

        default:
            Lo = gather_vpl(pr);      // MRCS, LightSlice, LightCuts
            break;
    }

    pr.val = Lo;
}

Rgb VirtualPointLight::connect_to_vpl(const Receiver &r, int pixel_index, int eye_bounce) {
    // receiver is not on light as we already check for that

    Rgb Lo;
    LocalGeometry dg(r);

    switch (estimate_mode) {
        case VplRadianceEstimate::LIGHTCUTS: 
        {
            Log::error() << "Lightcuts not supported." << endn;
            /*
            CutStat cut_stat;        
            Lo = lightcuts.cut(r, this, cut_stat);    
            Lo /= num_emitted_photons;
            */
            break;
        }

        case VplRadianceEstimate::MRCS: {
            Lo = mrcs.gather_clusters(r, this);
            Lo /= num_emitted_photons;
            break;
        }

        case VplRadianceEstimate::LIGHTSLICE: {
            Lo = lightslice.gather_clusters(r, this);
            Lo /= num_emitted_photons;
            break;
        }

        case VplRadianceEstimate::POWER_SAMPLING:
        case VplRadianceEstimate::NONE:
        default:
        {
            
            if (light_index == -1) return DefaultRgb::black;

            //if (all_vpls[light_index].get_bounce() != 0 || scene->get_integrator()->is_direct_lighting()) 
            //    return DefaultRgb::black;        

            BrdfPointLight *light = &all_vpls[light_index];
            if (light->get_bounce() + eye_bounce > max_bounce) return DefaultRgb::black;

            Lo = gather_point(scene, dg, r.wo, r.m, &all_vpls[light_index]);

            if (is_bias_compensation)
                Lo += kollig_compute_bias(r);
        
            if (estimate_mode == VplRadianceEstimate::NONE) {                
                Float progressive_correction_ratio = num_emitted_photons / total_bounce1;
                Float correction_factor = 1.0f / (num_evaluated_vpls * progressive_correction_ratio);

                Lo *= correction_factor;
            }
            if (estimate_mode == VplRadianceEstimate::POWER_SAMPLING) {
                Lo /= (pdf_light * num_emitted_photons);
            }

            break;
        }
    }

    return Lo;    
}

void VirtualPointLight::on_frame_end() {
    frame_index++;

    switch (estimate_mode) {
        case VplRadianceEstimate::BIAS_ONLY: {
            num_evaluated_bias_samples++;
            break;
        }

        case VplRadianceEstimate::NONE: {
            num_evaluated_vpls++;
            while (true) {
                light_index++;

                if (light_index >= all_vpls.size()) {
                    light_index = -1;
                    next_batch = true;
                    break;
                }                
            }
            break;
        }

        case VplRadianceEstimate::POWER_SAMPLING: {
            // choose a VPL for all receivers based on power distribution        
            Random &rd = *scene->get_random();
            powers.sample(rd, light_index, pdf_light);        

            Log::info() << "Light index: " << light_index << " " << pdf_light << endn;
            while (pdf_light <= 0.0f) {
                powers.sample(rd, light_index, pdf_light);
                Log::error() << "Resampling..." << endn;
            }
        
            num_power_samples_so_far++;
            if (num_power_samples_so_far > num_power_sampling_samples) {
                next_batch = true;
            }
            break;
        }

        case VplRadianceEstimate::MRCS:
        case VplRadianceEstimate::LIGHTCUTS: {
            next_batch = true;
            break;
        }
        
        case VplRadianceEstimate::LIGHTSLICE: 
        case VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE:
        {
            lightslice.on_frame_end();
            break;
        }
    }

    if (estimate_mode == VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE ||
        estimate_mode == VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE) {

        if (frame_index == 1 || (frame_index % 5 == 0)) {
            ostringstream oss;
            oss << scene->get_output_folder() << "/" << scene->get_name() << "_sample" << frame_index << ".exr";
            ImageFloat img(*scene->get_frame_buffer()->get_last_complete_frame_buffer());

            string file = oss.str();
            if (img.save(file.c_str()))
                Log::info() << file << " saved." << endn;
        }
    }

    if (next_batch) {
        // should reduce only when the next batch of VPL is generated (same as PPM)
        // relax clamp threshold for consistency
        iteration++;    
        clamp_threshold = base_clamp_threshold * pow(iteration, clamp_relax);
        Log::info() << "New clamping theshold: " << clamp_threshold << endn;
        next_batch = false;

        // now generate new VPLs
        // note that this function might evaluate some VPLs, so must set new clamping threshold before
        spawn_vpls();
    }
}

void VirtualPointLight::visualize_vpl(Scene *scene, ImageFloat *img, BrdfPointLight *light, Rgb color) {
    Camera *camera = scene->get_camera();
        
    BrdfPointLight *pl = light;

    Vec3 p = pl->org();
    Vec2 pixel;
    Vec3 wo_e;
    if (! camera->hit(p, pixel, wo_e)) return;

    img->accumulate(pixel, color);
}

/**
 *  Display the dA approximated proportionally to cosine(eye ray, normal at surface) as in Nayar's Shape from Interreflection paper.
 */
static void visualize_dA(Scene *scene, const Ray &r, Rgb &radiance) {
    Lights &lights = *scene->get_lights();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();

    HitRecord rec;
    
    if (! agg->hit(r, tmin, max_tmax, tick, rec)) return;

    if (rec.light) return;

    Vec3 p = r.org() + rec.t * r.dir();
    Vec3 pn = rec.normal;
    Material *m = rec.material;
    Vec3 wo = unit_vector(-r.dir());

    // display dA
	Float cosine = dot(-r.dir(), pn);
	if (cosine < 0) cosine = 0;
	radiance = Rgb(cosine, cosine, cosine);
	return;
}

// ----------------------------------------------------------------------------
// Kollig's compensation
// ----------------------------------------------------------------------------
Rgb VirtualPointLight::gather_point(Scene *scene, const LocalGeometry &dg, const Vec3 &wo, Material *m, const BrdfPointLight *vpl) {
    bool visible;
    return gather_point(scene, dg, wo, m, vpl, visible);
}

Rgb VirtualPointLight::gather_point(Scene *scene, const LocalGeometry &dg, const Vec3 &wo, Material *m, const BrdfPointLight *vpl, bool &visible) {    
    if (!m) return DefaultRgb::black;

    Vec3 p = dg.p;
    Vec3 pn = dg.n;
    Vec3 q = vpl->org();
    Vec3 qn = vpl->normal();        
    
    Aggregate *agg = scene->get_aggregate();    
    Float tmin     = scene->get_tmin();
    Float tick     = scene->get_tick();
    if (agg->hit(p, q, tmin, tick)) {
        visible = false;
        return DefaultRgb::black;
    }
    visible = true;

    Float G = FormFactor::form_factor_abs(p, pn, q, qn);    
    Vec3 wi = unit_vector(q - p);
    Rgb brdf = m->eval(dg, wo, wi);        
    Rgb light_brdf = vpl->brdf(-wi);    
    Rgb light_power = vpl->power();
    RadianceTerms term(light_power, brdf, light_brdf, G);
    
    clamp_term(term, clamp_threshold, clamp_mode);

    return term.radiance;     
}

/**
 * Find compensation for a gather point.
 *
 * TODO: For progressive, cache the gather point's path tracing, and perform per bounce for all gather points.
 */
Rgb VirtualPointLight::kollig_compute_bias(const Receiver &r) {
    if (! r.m) return DefaultRgb::black;

    Random &rd      = *scene->get_random();    
    Aggregate *agg  = scene->get_aggregate();
    Float tmin      = scene->get_tmin();
    Float max_tmax  = scene->get_max_tmax();
    Float tick      = scene->get_tick();

    LocalGeometry dg(r);
    ShadingGeometry sg(r);
    Vec3 p = dg.p;
    Vec3 pn = dg.n;
    
    Vec3 wi;
    Float pdf_wi;
    Rgb brdf = r.m->sample(rd, dg, r.wo, wi, pdf_wi);
    if (pdf_wi <= 0.0f) return DefaultRgb::black;

    // find the next intersection point    
    Ray ray(p, wi);
    HitRecord rec;
    if (! agg->hit(ray, tmin, max_tmax, tick, rec)) return DefaultRgb::black;

    if (rec.light) {
        return DefaultRgb::black;        
    }

    LocalGeometry dg_q(ray, rec);
    Receiver recv(ray, rec);
    Vec3 q = dg_q.p;
    Vec3 qn = dg_q.n;
    Vec3 unit_pq = unit_vector(q - p);
    Float G = FormFactor::form_factor_abs(p, pn, q, qn);
    
    // these will cancel with the form factor, but since we need form factor to examine clamping,
    // we need to save these terms in the power
    Float pdf_q = pdf_wi * fabs(dot(unit_pq, qn)) / (q - p).squared_length();
    Rgb contrib = 1.0f / pdf_q;
    RadianceTerms partial_term(contrib, brdf, DefaultRgb::white, G);

    if (clamp_mode == VplRadianceClamp::FORM_FACTOR) {
        RadianceTerms complete_term = partial_term;
        clamp_term(complete_term, clamp_threshold, clamp_mode);
        if (complete_term.is_clamped == false) {        // early termination
            return DefaultRgb::black;
        }
    }

    // continue path tracing from q to determine the BRDF at q and the power contribution to q
    // this completes the radiance term and allows us to know if clamping occurs.
    // this is general as we can clamp reflectivity or even entire illumination.
    Rgb bias = kollig_compute_bias_tracing(recv, partial_term, 2, max_bounce);
    return bias;
}

/**
 * This is very similar to trace_wellington() function but with additional checks
 * when connecting to light sources.
 */
Rgb VirtualPointLight::kollig_compute_bias_tracing(const Receiver &r, RadianceTerms &term, int bounce, int max_bounce) {
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

    if (bounce > max_bounce) return DefaultRgb::black;

    if (r.light) {
        // ignore double count of direct illumination
        return DefaultRgb::black;
    }

    Material *m = r.m;
    LocalGeometry dg(r);
    ShadingGeometry sg(r);
    
    Rgb bias;

    // connect to light and check clamping when the path is complete
    Float light_pdf;
    int k = scene->sample_light(rd(), light_pdf);
    if (k >= 0) {

        LightSamplingRecord sr;

        Rgb Lin = lights[k]->irradiance(scene, dg, sr) / light_pdf;
        if (sr.is_valid) {
        
            RadianceTerms complete_term = term;   

            if (bounce == 2) {
                complete_term.power *= Lin;     // the irradiance includes the cosine sg.cosine(-sr.wo)
                complete_term.light_brdf = sg.eval(-sr.wo);
            } else {
                complete_term.power *= Lin * sg.eval(-sr.wo);   
            }

            // TODO: clamp radiance must be evaluated on the path as if it is generated from the light source. 
            // because the probability of the path from eye or light is different, hence the contribution of a path 
            // generated by eye or light is different.
            // For now, we use reflectivity clamp, which is independent of how the path is generated.

            clamp_term(complete_term, clamp_threshold, clamp_mode);
            //bias = complete_term.radiance_no_clamp;
            if (complete_term.is_clamped) {
                
                // clamping only occurs at the second path segment from eye
                //bias = (1.0f - complete_term.w_clamp) * complete_term.radiance_no_clamp;                
                bias = (1.0f - complete_term.w_clamp);
            }

            if (bounce == 2 && clamp_mode == VplRadianceClamp::REFLECTIVITY &&         // early termination
                complete_term.is_clamped == false) {

                return DefaultRgb::black;
            }
        }
    }

    // continue tracing
    HitRecord rec;    
    Vec3 wi;
    Float pdf_wi;
    Rgb brdf = sg.sample(rd, wi, pdf_wi);
    if (pdf_wi == 0.0f || brdf == DefaultRgb::black) return bias;

    Ray ray(p, wi);
    if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
        if (rec.light)  // avoid double count
            return bias;

        // update partial terms
        if (bounce == 2) {
            term.power *= sg.cosine(wi) / pdf_wi;
            term.light_brdf = brdf;
        } else {
            term.power *= brdf * sg.cosine(wi) / pdf_wi;
        }

        Receiver recv(ray, rec);
        bias += kollig_compute_bias_tracing(recv, term, bounce + 1, max_bounce);
    }
    return bias;
}

Rgb VirtualPointLight::gather_bias(const PrimaryRay &pr) {
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick = scene->get_tick();
    Random &rd = *scene->get_random();
        
    Ray r = pr.r;
    HitRecord rec;
    if (!agg->hit(r, tmin, max_tmax, tick, rec)) return DefaultRgb::black;
        
    Receiver recv(r, rec);
    Rgb bias = kollig_compute_bias(recv);

    return bias;
}

Rgb VirtualPointLight::gather_vpl(const PrimaryRay &pr) {
        
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick = scene->get_tick();
    Random &rd = *scene->get_random();

    Rgb Lo = DefaultRgb::black;
    
    int bounce = 1;
    int max_bounce = this->get_max_bounce();

    Ray r = pr.r;
    HitRecord rec;
    if (!agg->hit(r, tmin, max_tmax, tick, rec)) return DefaultRgb::black;

    LocalGeometry dg(r, rec);
    Vec3 p = dg.p;
    Vec3 pn = dg.n;
    Vec3 wo = unit_vector(-r.dir());
        
    Receiver recv(r, rec);
        
    int width = scene->get_image_size().width;
    int pixel_index = pr.s.y * width + pr.s.x;
    
    // VPL does not store direct illumination
    bounce = 1;
    if (this->is_direct_lighting())             
        Lo += FirstBounce::radiance(scene, recv);

    Material *m = rec.material;
    if (m == NULL) return Lo;
                
    // just call gather cluster once to build the cache
    //connect_to_vpl(recv, pixel_index, bounce);
    Lo += connect_to_vpl(recv, pixel_index, bounce);
    return Lo;
}


Rgb VirtualPointLight::guided_path_trace(const PrimaryRay &pr) {
        
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick = scene->get_tick();
    Random &rd = *scene->get_random();

    Rgb Lo = DefaultRgb::black;
    
    int bounce = 1;
    int max_bounce = this->get_max_bounce();

    Ray r = pr.r;
    HitRecord rec;
    if (!agg->hit(r, tmin, max_tmax, tick, rec)) return DefaultRgb::black;

    LocalGeometry dg(r, rec);
    Vec3 p = dg.p;
    Vec3 pn = dg.n;
    Vec3 wo = unit_vector(-r.dir());
        
    Receiver recv(r, rec);
    
    return trace_wellington(scene, recv, 1, max_bounce);    
}

Rgb VirtualPointLight::trace_wellington(Scene *scene,
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
    
    HitRecord rec;    
    Vec3 wi;
    Float pdf_wi;
    
    Rgb radiance;
    Rgb brdf;

    if (mis_type == METRO_RANDOM) {
    
        bool is_brdf_sampling = rd() < 0.5f;    

        // direct lighting
        if (bounce > 1 || (bounce == 1 && this->is_direct_lighting())) {
            
            // TODO: here we have not yet perform MIS between BRDF and LightSlice technique for direct lighting
            // To implement this, we need to modify first_bounce() of the area light, 
            // so that when we can choose BRDF or LightSlice to sample a direction before combining with the light sample.
            //
            // This means BRDF and LightSlice are treated with weight 0.5 each, and multiply with 1 / pdf of BRDF/LightSlice selection,
            // which results in 1 in total.
            /*
            if (is_brdf_sampling) {
                radiance = FirstBounce::radiance(scene, r, NULL);
            } else {                
                DirectionSampler *dir_sampler = &lightslice;        // guided direct lighting
                radiance = FirstBounce::radiance(scene, r, dir_sampler);
            }*/

            radiance = FirstBounce::radiance(scene, r, NULL);
        }

        if (bounce >= max_bounce) return radiance;

        // indirect lighting
        if (is_brdf_sampling) {
            brdf = sg.sample(rd, wi, pdf_wi);            
        } else {
            brdf = lightslice.sample(rd, r, wi, pdf_wi);
        }
    
        if (pdf_wi > 0.0f) {

            Ray ray(p, wi);
            if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
                if (rec.light)  // avoid double count
                    return radiance;

                Receiver recv(ray, rec);
                Rgb Li = trace_wellington(scene, recv, bounce + 1, max_bounce);

                Float pdf_wi_other;
                if (is_brdf_sampling) {
                    //pdf_wi_other = lightslice.pdf(recv, wi);
                    // BUG found on Apr 07: must be r
                    pdf_wi_other = lightslice.pdf(r, wi);
                }
                else
                    pdf_wi_other = r.m->pdf(dg, r.wo, wi);

                Float weight = pdf_wi / (pdf_wi + pdf_wi_other);
                radiance += weight * Li * brdf * sg.cosine(wi) / pdf_wi * 2;
            }
        }
    
    } else if (mis_type == METRO_ONLY) {
        
        // debug patchy results
        //radiance = lightslice.get_slice_color(r);
        //radiance = lightslice.get_slice_sum(r);
        //return radiance;

        // test sampling
        //return lightslice.sample(rd, r, wi, pdf_wi);

        // use only Metropolis sampling   
        if (bounce > 1 || (bounce == 1 && this->is_direct_lighting())) {
            //DirectionSampler *dir_sampler = &lightslice;        // guided direct lighting
            //radiance = FirstBounce::radiance(scene, r, dir_sampler);

            radiance = FirstBounce::radiance(scene, r, NULL);
        }

        if (bounce >= max_bounce) return radiance;

        brdf = lightslice.sample(rd, r, wi, pdf_wi);
    
        if (pdf_wi > 0.0f) {

            Ray ray(p, wi);
            if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
                if (rec.light)  // avoid double count
                    return radiance;

                Receiver recv(ray, rec);
                Rgb Li = trace_wellington(scene, recv, bounce + 1, max_bounce);
                        
                radiance += Li * brdf * sg.cosine(wi) / pdf_wi;            
            }
        }
    }

    /***********************************
     * Integrate Vorba's LibImportance
     ***********************************/
    else if (mis_type == VORBA_ONLY) {

        if (bounce > 1 || (bounce == 1 && this->is_direct_lighting())) {
            radiance = FirstBounce::radiance(scene, r, NULL);
        }

        if (bounce >= max_bounce) return radiance;

        Vorba::PhotonSampler ps = importance_manager->get_sampler(LocalGeometry(r));
        brdf = ps.sample(rd, r, wi, pdf_wi); 

        if (pdf_wi > 0.0f) {

            Ray ray(p, wi);
            if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
                if (rec.light)  // avoid double count
                    return radiance;

                Receiver recv(ray, rec);
                Rgb Li = trace_wellington(scene, recv, bounce + 1, max_bounce);

                radiance += Li * brdf * sg.cosine(wi) / pdf_wi;
            }
        }

    }
    else if (mis_type == VORBA_RANDOM) {
        bool is_brdf_sampling = rd() < 0.5f;

        if (bounce > 1 || (bounce == 1 && this->is_direct_lighting())) {
            radiance = FirstBounce::radiance(scene, r, NULL);
        }

        if (bounce >= max_bounce) return radiance;

        Vorba::PhotonSampler ps = importance_manager->get_sampler(LocalGeometry(r));
        if (is_brdf_sampling) {
            brdf = sg.sample(rd, wi, pdf_wi);
        }
        else {

            brdf = ps.sample(rd, r, wi, pdf_wi);
        
        }

        if (pdf_wi > 0.0f) {

            Ray ray(p, wi);
            if (agg->hit(ray, tmin, max_tmax, tick, rec)) {
                if (rec.light)  // avoid double count
                    return radiance;

                Receiver recv(ray, rec);
                Rgb Li = trace_wellington(scene, recv, bounce + 1, max_bounce);

                Float pdf_wi_other;
                if (is_brdf_sampling) {
                    pdf_wi_other = ps.pdf(r, wi);
                }
                else
                    pdf_wi_other = r.m->pdf(dg, r.wo, wi);

                Float weight = pdf_wi / (pdf_wi + pdf_wi_other);
                radiance += weight * Li * brdf * sg.cosine(wi) / pdf_wi * 2;
            }
        }
    }

    return radiance;
}


Rgb VirtualPointLight::photon_gather(const Receiver &r, int pixel_index, int eye_bounce, int max_bounce) {
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
