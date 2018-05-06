#ifndef _IMAGE_VIEW_H_
#define _IMAGE_VIEW_H_

#include "math3.h"
#include "gl_view.h"
#include "boundingbox.h"
#include "integrator.h"
#include "stats.h"
#include "frame.h"
#include "events.h"

#include <boost/thread.hpp>
namespace Renzoku {

struct ByteBlock;

/**
 * A simple block-based image viewer that draws the completed blocks 
 * by MonteCarloIntegrator onto screen.
 */
class ImageBlockView : public GLView, public FrameObserver {
public:
    ImageBlockView(Scene *scene);
    ImageBlockView(Scene *scene, Size2 window_size);
    ~ImageBlockView();

    void init();
    void reset();

    void reshape(int width, int height);
    void display();
    void keyboard(unsigned char key, int x, int y);
    void mouse(MouseButton button, MouseState state, int x, int y);
    void motion(int x, int y);

    void on_show();
    void on_hide();
    void on_close();
    
    void save_intermediate(bool save);

    virtual void on_complete_frame(int tid, ImageFloat *img);
    virtual void on_update_image(int tid, ImageByte *img, int x0, int y0, int x1, int y1);
    virtual void on_update_border(int tid, int x0, int y0, int x1, int y1);

private:
    void init(Scene *scene, Size2 window_size);

    void save_linear_last_complete_frame(const char *file);    
    void save_linear_screen(const char *file);

    inline ByteBlock * block_queue_push();
    inline ByteBlock*  block_queue_front();
    inline void        block_queue_pop();

private:
    FrameBuffer *frame;
    ImageByte *byte_buffer;

    /**
     * Queue of blocks to render to screen
     */
    ByteBlock **block_queue;
    int num_blocks_in_queue;
    int max_blocks_in_queue;
    int queue_first, queue_last;

    int win_width;
    int win_height;

    string tmp_file_template;           // temporary file name
    string tmp_screen_template;         // temporary file name for screenshot
    bool is_save_intermediate;

    bool redraw_all;
    ByteBlock *display_block;

    Stats fps_counter;
    int cur_frame;

    int num_colors;
    Byte reds[32];
    Byte greens[32];
    Byte blues[32];

private:
    //boost::mutex mtx;

public:
    MouseObservable mouse_event;
};

struct ByteBlock {
    Byte *data;
    int capacity;
    int x0, y0, x1, y1;

    ByteBlock(int capacity) {
        data = new Byte[capacity];
        x0 = y0 = x1 = y1 = 0;
        this->capacity = capacity;
    }

    void set_bound(int x0, int y0, int x1, int y1) {
        this->x0 = x0;
        this->y0 = y0;
        this->x1 = x1;
        this->y1 = y1;
    }

    void copy(ByteBlock *src) {        
        memcpy(data, src->data, src->capacity * sizeof(Byte));
        x0 = src->x0;
        y0 = src->y0;
        x1 = src->x1;
        y1 = src->y1;
        capacity = src->capacity;
    }

    void reset() {
        memset(data, 0, sizeof(Byte) * capacity);
    }
};

inline ByteBlock* ImageBlockView::block_queue_push() {
    if (num_blocks_in_queue < max_blocks_in_queue) {
        ByteBlock *last  = block_queue[queue_last];
        num_blocks_in_queue++;
        queue_last = (queue_last + 1) % max_blocks_in_queue;
        return last;
    }
    return NULL;
}

inline ByteBlock* ImageBlockView::block_queue_front() {
    if (num_blocks_in_queue > 0) {
        return block_queue[queue_first];
    }
    return NULL;
}

inline void ImageBlockView::block_queue_pop() {
    if (num_blocks_in_queue > 0) {
        num_blocks_in_queue--;
        queue_first = (queue_first + 1) % max_blocks_in_queue;
    }
}

} // end namespace

#endif