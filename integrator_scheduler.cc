#include "integrator_scheduler.h"
#include "scene.h"
#include "image.h"
#include "camera.h"
#include "frame.h"
#include "ray.h"

#include "monte_carlo_integrator.h"
#include "stats.h"

#include <boost/thread.hpp>

namespace Renzoku {
        
static boost::mutex scheduler_mtx;
static boost::thread_group scheduler_thread_group;

IntegratorScheduler::IntegratorScheduler() {
}

void IntegratorScheduler::initialize(Scene *scene) {
    is_threads_launched = false;
    num_threads = 1;

    this->scene = scene;
    this->frame = scene->get_frame_buffer();
    this->integrator = dynamic_cast<MonteCarloIntegrator *>(scene->get_integrator());
    if (this->integrator == NULL) {
        throw BaseException("Monte Carlo integrator expected.");
    }
}

void IntegratorScheduler::sample_pack_to_primary_ray_bundle(const SensorSamplePacks &packs, PrimaryRayBundles &bundles) {
    Camera *camera = scene->get_camera();
    for (int i = 0; i < packs.size(); ++i) {
        SensorSamplePack pack = packs[i];
        
        PrimaryRayBundle bundle;
        bundle.x0 = pack.x0;        
        bundle.y0 = pack.y0;
        bundle.x1 = pack.x1;
        bundle.y1 = pack.y1;
        for (int j = 0; j < pack.samples.size(); ++j) {
            SensorSample s = pack.samples[j];
            Ray r = camera->shoot_ray(s);
            PrimaryRay pr(r, s);
            bundle.rays.push_back(pr);
        }
        bundle.tag = (PrimaryRayBundle::RayBundleEnum)pack.tag;
        bundles.push_back(bundle);
    }
}

void IntegratorScheduler::start() {
    SensorSamplePacks sample_packs;
    PrimaryRayBundles bundles;
    Random &rd = *scene->get_random();
        
    Log::info() << "Threads: " << num_threads << endn;
    
    sample_packs.clear();
    Sensor *sensor = scene->get_camera()->get_sensor();
    sensor->reset(); 
    
    if (sensor->has_debug_pixels()) {
        sensor->get_debug_sample_packs(sample_packs);
        sample_pack_to_primary_ray_bundle(sample_packs, bundles); 
        this->sequential(&bundles);
        return;
    }
    
    const int num_packs_per_round = num_threads;
    int cur_samples = 0;
    Stats stats;
    float frame_time = 0.f;
    stats.tic();
    while (cur_samples < integrator->get_num_samples()) {
        if (bundle_pool.size() > 128) {
            boost::this_thread::sleep(boost::posix_time::microseconds(1));
            continue;
        }

        bool has_end_frame;
        sample_packs.clear();
        sensor->get_sample_packs(num_packs_per_round, rd, sample_packs, has_end_frame);

        // generate primary rays by CPU
        bundles.clear();        
        sample_pack_to_primary_ray_bundle(sample_packs, bundles); 

        if (num_threads > 1) {
            this->queue(&bundles);
        } else {
            this->sequential(&bundles);
        }

        if (has_end_frame) {
            cur_samples++;
            stats.toc();
            frame_time += stats.elapsed();

            Log::info() << "Samples: " << cur_samples << "\t " << cur_samples * 1000 / frame_time << " samples/sec" << endn;            

            stats.tic();
        }
    }

    scheduler_thread_group.join_all();
    Log::info() << "Integration complete." << endn;
}

void IntegratorScheduler::on_primary_ray_bundle_start(int tid, PrimaryRayBundle &b) {
    integrator->observable.notify_update_border(tid, b.x0, b.y0, b.x1, b.y1);

    if (b.tag == PrimaryRayBundle::BUNDLE_BEGIN_FRAME) {
        integrator->on_frame_begin();
    }
}

void IntegratorScheduler::on_primary_ray_bundle_done(int tid, PrimaryRayBundle &b) {
    ImageFloat *img = scene->get_frame_buffer()->get_current_buffer();
    integrator->observable.notify_update_image(tid, img, b.x0, b.y0, b.x1, b.y1);

    if (b.tag == PrimaryRayBundle::BUNDLE_END_FRAME) {
        integrator->observable.notify_complete_frame(tid, img);

        // FIXME: keep track of total blocks so far per frame to make sure frame is complete in case of multithread.

        integrator->on_frame_end();
    }
}

void IntegratorScheduler::sequential(PrimaryRayBundles *bundles) {
    for (int i = 0; i < bundles->size(); ++i) {
        PrimaryRayBundle &pri = (*bundles)[i];
        
        this->on_primary_ray_bundle_start(0, pri);
        integrator->integrate(pri, frame);
        this->on_primary_ray_bundle_done(0, pri);
    }
}

void IntegratorScheduler::queue(PrimaryRayBundles *bundles) {
    for (int i = 0; i < bundles->size(); ++i) {
        PrimaryRayBundle &pri = (*bundles)[i];
        bundle_pool.push(pri);
    }

    // launch persistent threads
    if (! is_threads_launched) {
        is_threads_launched = true;
        if (scheduler_thread_group.size() > 0) {
            scheduler_thread_group.join_all();
        }
        
        for (int i = 0; i < num_threads; ++i) 
            hungry[i] = 0;

        for (int i = 0; i < num_threads; ++i) {
            boost::thread *t = new boost::thread(boost::bind(&IntegratorScheduler::worker, this, i));
            scheduler_thread_group.add_thread(t);
        }
    }   
}

    /*
    // each thread per ray and relaunch. Very slow.
    const int num_threads = 4;
    boost::thread **t = new boost::thread*[num_threads];
    for (int i = 0; i < b.rays.size(); i += num_threads) {
        for (int j = 0; j < num_threads; ++j) {
            int idx = i + j;
            t[j] = new boost::thread(call_radiance, this, b.rays[idx], boost::ref(Lo[idx]));
        }

        for (int j = 0; j < num_threads; ++j) {
            t[j]->join();

            delete t[j];
        }        
    }
    delete [] t;    
    */

void IntegratorScheduler::worker(int tid) {
    Stats stats;
    while (true) {
        //Log::info() << "Thread : " << tid << endn;
        
        stats.tic();
        if (bundle_pool.empty()) {
            boost::this_thread::sleep(boost::posix_time::microseconds(1));
            //Log::debug() << "-- Thread " << tid << " starved." << endn;

            // if starved a few times then stop
            hungry[tid]++;
            if (hungry[tid] > MAX_HUNGRY) {
                //Log::debug() << "-- Thread " << tid << " ended." << endn;
                break;
            }

            continue;
        }
        stats.toc();
        //Log::info() << "-- " << tid << " -- Pool check : " << stats.elapsed() << endn;

        stats.tic();
        scheduler_mtx.lock();
        PrimaryRayBundle pri = bundle_pool.front();
        bundle_pool.pop();
        scheduler_mtx.unlock();
        stats.toc();
        //Log::info() << "-- " << tid << " -- Get bundle : " << stats.elapsed() << endn;

        stats.tic();
        this->on_primary_ray_bundle_start(tid, pri);
        stats.toc();

        stats.tic();
        integrator->integrate(pri, frame);
        stats.toc(); 
        //Log::info() << "-- " << tid << " -- Ray trace : " << stats.elapsed() << endn;

        stats.tic();
        this->on_primary_ray_bundle_done(tid, pri);
        stats.toc();
        //Log::info() << "-- " << tid << " -- Notification : " << stats.elapsed() << endn;
    }
}

} // end namespace
