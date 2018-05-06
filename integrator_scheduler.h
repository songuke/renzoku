#ifndef _MONTE_CARLO_INTEGRATOR_SCHEDULER_H_
#define _MONTE_CARLO_INTEGRATOR_SCHEDULER_H_

#include "common.h"
#include "size.h"
#include "image.h"
#include "observable.h"
#include "sensor.h"

#include <queue>
using namespace std;

namespace Renzoku {

const int MAX_THREADS = 16;
const int MAX_HUNGRY = 128;

class IParallelIntegrator {
public:
    inline int get_num_threads() const;
    inline void set_num_threads(int threads);

protected:
    int num_threads;
};

inline int IParallelIntegrator::get_num_threads() const {
    return num_threads;
}

inline void IParallelIntegrator::set_num_threads(int threads) {
    num_threads = threads;
    if (num_threads > MAX_THREADS) {
        num_threads = MAX_THREADS;
        Log::info() << "Max number of threads supported is " << MAX_THREADS << endn;
    }
}

class MonteCarloIntegrator;

/**
 * The scheduler determines which parts on the image plane should be integrated
 * and schedule for the integration. 
 * 
 * The integrator should not know about the rendering schedule or the display process.
 */
class IntegratorScheduler : public IParallelIntegrator {

public:
    IntegratorScheduler();
    
    void initialize(Scene *scene);
    void start();

protected:
    /**
     * Simply unpack each ray and trace it. No parallel.
     */
    void sequential(PrimaryRayBundles *bundles);

    /**
     * Queue each ray bundle to a pool for each worker thread to process.
     */
    void queue(PrimaryRayBundles *bundles);

private:
    void on_primary_ray_bundle_start(int tid, PrimaryRayBundle &bundle);
    void on_primary_ray_bundle_done(int tid, PrimaryRayBundle &bundle);
    void sample_pack_to_primary_ray_bundle(const SensorSamplePacks &packs, PrimaryRayBundles &bundles);

private:
    void worker(int thread_index);

private:
    Scene *scene;
    FrameBuffer *frame;
    MonteCarloIntegrator *integrator;

    bool is_threads_launched;
    std::queue<PrimaryRayBundle> bundle_pool;       // keep copy of ray bundles

    int hungry[MAX_THREADS];        
};

} // end namespace

#endif