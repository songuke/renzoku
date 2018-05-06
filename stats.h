#ifndef _STATS_H_
#define _STATS_H_

#include <ctime>
#include <queue>
using std::clock_t;

namespace Renzoku {

class Stats {
public:
    Stats();

    int tic();
    int toc();

    /**
     * Return the time elapsed (in ms) between last toc and tic.
     */
    float elapsed() const;

    /**
     * Report estimated frame rate over the history of time elapsed.
     */ 
    float fps();

    float get_total_elapsed();

    inline int get_total_frames() const;

    void reset_total_elapsed();

    void print_elapsed_milliseconds();
    void print_elapsed_seconds();

protected:
    clock_t last_tic;
    clock_t last_toc;

    /**
     * We estimate the frame rate using the following formula:
     * FPS = n / T(n)           (total n frames divided by the running time of n frames
     *     ~ a * (n - 1) / T(n - 1) + (1 - a) * 1 / t(n)
     *     ~ a * old_FPS + (1 - a) * immediate_FPS
     *     ~ a * old_FPS * (1 - a) * m / t(m)      (m is number of frames between each FPS estimation call.)
     *
     * where T(n) = sum(t(i))    i = 1 .. n.
     *
     * The weight a determines how much the previous frame rate estimation affects the frame rate estimation. 
     */
    float last_fps;
    const float a;
    const float one_minus_a;

    /**
     * An accumulated timer
     */
    clock_t total_clocks;
    int total_frames;
};

inline int Stats::get_total_frames() const {
    return total_frames;
}

}  // end namespace

#endif
