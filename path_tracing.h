#ifndef _PATH_TRACING_H_
#define _PATH_TRACING_H_

#include "common.h"
#include "monte_carlo_integrator.h"

namespace Renzoku {

class PathTracing : public MonteCarloIntegrator {
public:
    PathTracing();

    virtual Rgb radiance(const Receiver &);
    virtual void initialize(Scene *scene);

    /**
     * Export the trace function for VPL sampling use.
     */
    static Rgb trace(Scene *scene, const Receiver &r, int max_bounce);

    virtual void on_frame_end();

protected:
    /**
     * Experimental
     */
    void radiance_field(const Receiver &r);
};

} // end namespace Renzoku

#endif
