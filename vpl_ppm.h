#ifndef _VIRTUAL_POINT_LIGHT_PPM_H_
#define _VIRTUAL_POINT_LIGHT_PPM_H_

#include "monte_carlo_integrator.h"
#include "kdtree.h"
#include "vpl_interface.h"
#include "vpl_generator.h"

namespace Renzoku {

struct DomainStats {
    Rgb acc_radiance;
    long long num_photons;
    Float radius;

    DomainStats() 
        : num_photons(0), acc_radiance(DefaultRgb::black), radius(-1.0f) // negative
    {

    }
};

class VplPpm : public MonteCarloIntegrator, 
               public IntegratorObserver,
               public IVirtualPointLight {
public:
    VplPpm();
    ~VplPpm();
    virtual void initialize(Scene *scene);

    /**
     * For multi-pass photon mapping.
     */
    inline void set_nearest_photons(int num);

    /**
     * For progressive photon mapping.
     */
    inline void use_progressive_photon(bool sppm);    

    /**
     * In order to resolve more details, number of photons must be increased, 
     * and radius must be decreased.
     * 
     * alpha determines the radius reduction rate.
     *
     * alpha = 1 yields no radius reduction, and thus is equivalent to multi-pass photon mapping.
     */
    inline void set_alpha(Float alpha);
        
    virtual void radiance(PrimaryRay &r);
    
    virtual void on_update_image(int tid, ImageFloat *img, int x0, int y0, int x1, int y1) {}
    virtual void on_complete_frame(int tid, ImageFloat *img);

protected:
    /**
     * Multi-pass photon gathering (traditional photon maping).
     *
     * This only computes the average radiance of multi-pass using
     * k-nearest photons.
     *
     * This approach cannot resolve details (e.g., caustics). It only smooths out
     * the artifacts; the caustics might still be blurry.
     */
    virtual Rgb photon_gather(const Receiver &r, int eye_bounce, int max_bounce);

    /**
     * Stochastic progressive photon mapping.
     *
     * This approach keeps track of photon statistics at each pixel.
     *
     * Photons of which bounce is greater than max bounce is ignored. 
     * (for final gathering use)
     */
    virtual Rgb photon_gather(const Receiver &r, int pixel_index, int eye_bounce, int max_bounce);

protected:

    /**
     * Gather at the first non-specular vertex.
     *
     * This ensures mirror-like surface at the first hit point to be rendered correctly.
     * (Photon gather at mirror-like surface is not efficient.)
     *
     */
    Rgb mixed_gather(const PrimaryRay &pr);

    
protected:
    KdTree *kdtree;
    int num_nearest_photons;

    bool use_sppm;
    DomainStats *domain_stats;    
    long long num_emitted_photons;
    Float alpha;
};

inline void VplPpm::set_nearest_photons(int num) {
    num_nearest_photons = num;
}

inline void VplPpm::use_progressive_photon(bool sppm) {
    use_sppm = sppm;
}

inline void VplPpm::set_alpha(Float alpha) {
    this->alpha = alpha;
}

} // end namespace Renzoku

#endif

