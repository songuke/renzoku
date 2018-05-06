#ifndef _VIRTUAL_POINT_LIGHT_H_
#define _VIRTUAL_POINT_LIGHT_H_

#include "common.h"
#include "monte_carlo_integrator.h"
#include "matrix.h"
#include "camera.h"
#include "vpl_generator.h"
#include "vpl_interface.h"

//#include "lightcuts.h"
#include "mrcs.h"
#include "lightslice.h"

#include "other/vorba/photon_sampler.h"

#include "kdtree.h"
#include "vpl_ppm.h"

#include "events.h"

namespace Renzoku {

class VirtualPointLight : public MonteCarloIntegrator, 
                          public IVirtualPointLight,
                          public IVirtualPointLightEvaluator {
public:
    VirtualPointLight();
    ~VirtualPointLight();

    virtual void initialize(Scene *scene);

    virtual void radiance(PrimaryRay &pr);

    virtual Rgb radiance(Scene *scene, const Receiver &r, const BrdfPointLight &light, bool &visible);
    virtual Rgb incident_radiance(Scene *scene, const Receiver &r, const BrdfPointLight &vpl);
        
    virtual void on_frame_end();

    inline void set_clamping(VplRadianceClamp::Type mode);
    inline void set_clamping_threshold(Float clamp_threshold);
    inline void set_clamping_relaxation(Float alpha);
    
    void set_bias_compensation(bool bias);
    void set_bias_compensation_samples(int num_samples);

    void set_radiance_estimate(VplRadianceEstimate::Type mode);

    void set_power_sampling_samples(int samples);

protected:    
    void spawn_vpls();

    /**
     * Compute the bias energy due to clamping by following Kollig and Keller, 2004 paper.
     */
    Rgb  compute_bias_energy(Scene *scene,
                             const LocalGeometry &dg,
                             const Vec3 &wo, Material *m) const;

    
    /** 
     * Light sources are added in the last pass.
     */
    void render_area_light_sources(Scene *scene, const Ray &r, Rgb &radiance) const;
    
protected:
    /**
     * Write 3D coordinates of a set of VPLs to a .PLY file 
     */
    static void debug_vpl(BrdfPointLights &vpls, const char *file);
    
public:
    /**
     * Project VPL onto the result image.
     */
    static void visualize_vpl(Scene *scene, ImageFloat *img, BrdfPointLight *light, Rgb color);

protected:
    /**
     * Evaluate a VPL to a point
     */
    Rgb gather_point(Scene *scene, const LocalGeometry &dg, const Vec3 &wo, Material *m, const BrdfPointLight *vpl);
    Rgb gather_point(Scene *scene, const LocalGeometry &dg, const Vec3 &wo, Material *m, const BrdfPointLight *vpl, bool &visible);
        
    /**
     * Sample a new ray at the gather point. Do recursive compensation if the next point is not counted in VPL clamped illumination.
     */
    Rgb gather_bias(const PrimaryRay &pr);
    Rgb  kollig_compute_bias(const Receiver &r);
    Rgb  kollig_compute_bias_tracing(const Receiver &r, RadianceTerms &partial_term, int bounce, int max_bounce);
    
protected:
    VplRadianceClamp::Type clamp_mode;      /// perform clamping during gathering to discard bright spots. Default: true.
    Float clamp_threshold;
    Float base_clamp_threshold;
    Float clamp_relax;

    bool is_bias_compensation;              /// calculate and add back missing energy due to clamping. Default: false.
    int num_bias_compensation_samples;
    int num_evaluated_bias_samples;

    ImageFloat *img_vpl;                    /// result of evaluating clamped VPL
    ImageFloat *img_bias;                   /// bias energy image

    BrdfPointLights all_vpls;    
    int standard_offset;
    int num_evaluated_vpls;

    int frame_index;

    long long num_emitted_photons;
    long long total_bounce1;

    int pass;                               /// rendering pass. 1: VPL gathering. 2: bias compensation.
    
    int light_index;                        /// to replace standard_offset
    Float pdf_light;

    int iteration;
    bool next_batch;                        /// indicate spawning new VPLs
        
    VplRadianceEstimate::Type estimate_mode;

    DiscretePdf powers;                     /// power sampling
    int num_power_samples_so_far;
    int num_power_sampling_samples;         
    
public:
    //Lightcuts lightcuts;
    MatrixRowColumnSampling mrcs;

    LightSlice lightslice;    
protected:

    enum Mis {
        METRO_ONLY,
        METRO_RANDOM,
        VORBA_ONLY,
        VORBA_RANDOM
    } mis_type;
    
protected:
    Float incident_radiance_clamp_threshold;
    int num_virtual_indirect_samples;

public:
    inline void set_incident_radiance_clamping_threshold(Float threshold);
    inline void set_virtual_indirect_samples(int num);

protected:
    Vorba::ImportanceManager *importance_manager;

protected:
    Rgb gather_vpl(const PrimaryRay &pr);    

    Rgb guided_path_trace(const PrimaryRay &pr);    
    Rgb trace_wellington(Scene *scene,
                         const Receiver &r, int bounce, int max_bounce);

    Rgb connect_to_vpl(const Receiver &r, int pixel_index, int eye_bounce);

/*
 * Photon mapping port
 */
protected:
    Rgb photon_gather(const Receiver &r, int pixel_index, int eye_bounce, int max_bounce);

protected:
    KdTree *kdtree;
    int num_nearest_photons;

    DomainStats *domain_stats;    
    Float alpha;

};

inline void VirtualPointLight::set_clamping(VplRadianceClamp::Type mode) {
    this->clamp_mode = mode;
}

inline void VirtualPointLight::set_clamping_threshold(Float threshold) {
    clamp_threshold = threshold;
    base_clamp_threshold = threshold;
}

inline void VirtualPointLight::set_incident_radiance_clamping_threshold(Float threshold) {
    incident_radiance_clamp_threshold = threshold;
}

inline void VirtualPointLight::set_virtual_indirect_samples(int samples) {
    num_virtual_indirect_samples = samples;
}

inline void VirtualPointLight::set_clamping_relaxation(Float alpha) {
    clamp_relax = alpha;
}

inline void VirtualPointLight::set_bias_compensation(bool bias) {
    this->is_bias_compensation = bias;
}

inline void VirtualPointLight::set_bias_compensation_samples(int num_samples) {
    num_bias_compensation_samples = num_samples;
}

inline void VirtualPointLight::set_radiance_estimate(VplRadianceEstimate::Type mode) {
    estimate_mode = mode;

    switch (estimate_mode) {
    case VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE:
        //mis_type = METRO_ONLY;
        mis_type = METRO_RANDOM;
        break;
    case VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE:
        //mis_type = VORBA_ONLY; 
        mis_type = VORBA_RANDOM;
        break;
    }
}

inline void VirtualPointLight::set_power_sampling_samples(int samples) {
    num_power_sampling_samples = samples;
}

} // end namespace Renzoku

#endif

