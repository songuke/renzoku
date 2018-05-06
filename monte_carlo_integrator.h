#ifndef _MONTE_CARLO_INTEGRATOR_H_
#define _MONTE_CARLO_INTEGRATOR_H_

#include "common.h"
#include "math3.h"
#include "integrator.h"
#include "sensor.h"

namespace Renzoku {
       
class IMonteCarloIntegrator {
public:
    inline int get_num_samples() const;
    inline void set_num_samples(int num_samples);

protected:
    int num_samples;
};

inline int IMonteCarloIntegrator::get_num_samples() const {
    return num_samples;
}

inline void IMonteCarloIntegrator::set_num_samples(int num_samples) {
    this->num_samples = num_samples;
}

class MonteCarloIntegrator : public Integrator,
                             public IMonteCarloIntegrator {
public:
    MonteCarloIntegrator();

    virtual void initialize(Scene *scene); 
    virtual inline string get_suffix() const;

    /**
     * Perform integration of the bundle and provide progressive update to the frame buffer.
     */
    virtual void integrate(PrimaryRayBundle &b, FrameBuffer *frame);

    virtual void on_frame_begin() {}
    virtual void on_frame_end() {}

protected:
    virtual void radiance(PrimaryRay &r);
    virtual Rgb radiance(const Ray &r);
    virtual Rgb radiance(const Receiver &r);
    
protected:
    string suffix;
    
    ImageFloat *target;
    ImageFloat *weight;
};

inline string MonteCarloIntegrator::get_suffix() const {
    return suffix;
}

} // end namespace

#endif