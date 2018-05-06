#ifndef _FRAME_H_
#define _FRAME_H_

#include "image.h"
#include "tonemap.h"
#include "integrator.h"
#include "observable.h"

namespace Renzoku {

class FrameObserver {
public:    
    virtual void on_update_image(int tid, ImageByte *img, int x0, int y0, int x1, int y1) = 0;
    virtual void on_complete_frame(int tid, ImageFloat *img) = 0;

    /** 
     * For progress update.
     */
    virtual void on_update_border(int tid, int x0, int y0, int x1, int y1) {}
};
typedef vector<FrameObserver *> FrameObservers;


class FrameObservable : public Observable<FrameObserver> {
public:
    void notify_update_image(int tid, ImageByte *img, int x0, int y0, int x1, int y1) {
        FrameObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_update_image(tid, img, x0, y0, x1, y1);
    }

    void notify_complete_frame(int tid, ImageFloat *img) {
        FrameObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_complete_frame(tid, img);
    }

    void notify_update_border(int tid, int x0, int y0, int x1, int y1) {
        FrameObservers::iterator i;
        for (i = observers.begin(); i != observers.end(); ++i) 
            (*i)->on_update_border(tid, x0, y0, x1, y1);
    }
};

/**
 * Frame buffer design:
 *
 * Frame stores the latest progressive result, and the its tonemapped buffer. It also maintains the last complete frame buffer. 
 * Frame does not maintain accumulation and weight buffer for Monte Carlo estimation as they are specific to Monte Carlo methods.
 * Those are maintained by MonteCarloIntegrator class.
 *
 * Each frame is considered complete when a round of certain operations are finished for all pixels in the image (e.g., each pixel 
 * complete accumulating a sample.)
 */
class FrameBuffer : public IntegratorObserver {
public:
    FrameBuffer(Size2 img_size);
    ~FrameBuffer();

    void reset();

    void set_display_correction_factor(Float d) {
        display_correction_factor = d;
    }
    
    ImageFloat *get_current_buffer() const {
        return buffer;
    }
    
    ImageFloat *get_tonemap_buffer() const {
        return tonemap_buffer;
    }

    ImageByte *get_byte_buffer() const {
        return byte_buffer;
    }

    ImageFloat *get_last_complete_frame_buffer() const {
        return frame_buffer;
    }

    int get_index() const {
        return index;
    }

    virtual void on_update_image(int tid, ImageFloat *img, int x0, int y0, int x1, int y1);
    virtual void on_update_border(int tid, int x0, int y0, int x1, int y1);
    virtual void on_complete_frame(int tid, ImageFloat *img);    

private:    
    ImageFloat *buffer;                     // the latest result (progressive)
    ImageFloat *tonemap_buffer;             // tone mapped buffer to display the progressive result. 
    ImageFloat *frame_buffer;               // last complete frame buffer    
    ImageByte  *byte_buffer;                // final frame buffer. No processing should be done on this buffer to avoid display flickering.
    
    Float scale;                            // to linear scale tonemapped data to [0, 1]

    ReinhardToneMap *tonemapper;
    GammaToneMap *gamma_mapper;

    Float display_correction_factor;        // Default: 1.
    int        index;    

public:
    FrameObservable observable;
};

} // end namespace

#endif