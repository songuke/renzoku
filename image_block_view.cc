#include <GL/freeglut.h>

#include "image_block_view.h"
#include "scene.h"
#include "viewer.h"

namespace Renzoku {

ImageBlockView::ImageBlockView(Scene *scene) : GLView(scene), cur_frame(0) {
    Size2 img_size = scene->get_image_size();
    init(scene, img_size);
}

ImageBlockView::ImageBlockView(Scene *scene, Size2 window_size) : GLView(scene), cur_frame(0) {
    scene->set_image_view(this);
    init(scene, window_size);
}

void ImageBlockView::init(Scene *scene, Size2 window_size) {
    cur_frame = 0;
    
    frame          = scene->get_frame_buffer();
    byte_buffer    = frame->get_byte_buffer();     

    this->win_width  = window_size.width;
    this->win_height = window_size.height;

    tmp_file_template   = scene->get_name() + "_%d.exr";
    tmp_screen_template = scene->get_name() + "_snapshot_%d.exr";

    redraw_all = false;        
    display_block = new ByteBlock(3 * (int)window_size.height * (int)window_size.width);
    
    max_blocks_in_queue = 64;
    block_queue = new ByteBlock*[max_blocks_in_queue];
    for (int i = 0; i < max_blocks_in_queue; ++i)
        block_queue[i] = new ByteBlock(3 * (int)window_size.height * (int)window_size.width);
    num_blocks_in_queue = 0;
    queue_first = 0;
    queue_last = 0;

    // predefined colors
    num_colors = 8;
    for (int i = 0; i < num_colors; ++i) {
        Rgb c = Rgb::from_hsv((Float)(360 * i) / num_colors, 0.75f, 0.75f);
        reds[i] = c.r * 255;
        greens[i] = c.g * 255;
        blues[i] = c.b * 255;
    }
}

void ImageBlockView::init() {
}

void ImageBlockView::reset() {

}

void ImageBlockView::reshape(int width, int height) {
    this->win_width = width;
    this->win_height = height;

    redraw_all = true;
}

void ImageBlockView::on_show() {
    glViewport(0, 0, win_width, win_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0f, win_width, 0.0f, win_height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glDrawBuffer(GL_FRONT_AND_BACK);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    
    glFinish();             // must clear the screen as the first block might not arrive 
                            // (say in debug mode it can be slow.)

    frame->reset();    
    display_block->reset();

    // relaunch    
    scene->start();
}

void ImageBlockView::on_hide() {
}

void ImageBlockView::on_update_image(int tid, ImageByte *byte_buffer, int x0, int y0, int x1, int y1) {
    ByteBlock *back_block = block_queue_push();
    if (! back_block) return;
    
    int k = 0;
    for (int i = y0; i < y1; ++i) {
        for (int j = x0; j < x1; ++j) {
            back_block->data[k++] = (*byte_buffer)(i, j, 0);
            back_block->data[k++] = (*byte_buffer)(i, j, 1);
            back_block->data[k++] = (*byte_buffer)(i, j, 2);
        }
    }
    back_block->set_bound(x0, y0, x1, y1);

    ////mtx.lock();
    //mtx.unlock();
}

void ImageBlockView::on_update_border(int tid, int x0, int y0, int x1, int y1) {
    int k = 0;
    // copy current display data and put a border on it
    ByteBlock *back_block = block_queue_push();
    if (! back_block) return;
    
    for (int i = y0; i < y1; ++i) {
        for (int j = x0; j < x1; ++j) {
            back_block->data[k++] = (*byte_buffer)(i, j, 0);
            back_block->data[k++] = (*byte_buffer)(i, j, 1);
            back_block->data[k++] = (*byte_buffer)(i, j, 2);
        }
    }
    Byte red, green, blue;
    red     = reds[tid];
    green   = greens[tid];
    blue    = blues[tid];

    k = 0;
    int stride = 3 * (x1 - x0);
    for (int j = x0; j < x1; ++j) {                
        back_block->data[k++] = red;
        back_block->data[k++] = green;
        back_block->data[k++] = blue;
    }
    k = (y1 - y0 - 1) * stride;
    for (int j = x0; j < x1; ++j) {                
        back_block->data[k++] = red;
        back_block->data[k++] = green;
        back_block->data[k++] = blue;
    }
    k = 0;
    for (int i = y0; i < y1; ++i) {
        back_block->data[k] = red;
        back_block->data[k+1] = green;
        back_block->data[k+2] = blue;
        k += stride;
    }
    k = 3 * (x1 - x0 - 1);
    for (int i = y0; i < y1; ++i) {
        back_block->data[k] = red;
        back_block->data[k+1] = green;
        back_block->data[k+2] = blue;
        k += stride;
    }
    back_block->set_bound(x0, y0, x1, y1);
}

void ImageBlockView::on_complete_frame(int tid, ImageFloat *img) {

}

void ImageBlockView::display() {
    // transfer display buffer image to byte buffer    
    Byte *byte_data = byte_buffer->get_data();    
    int img_width = byte_buffer->get_width();
    int img_height = byte_buffer->get_height();

    // model view matrix must be reset each round to position blocks properly
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // disable depth test so blocks can be overwritten at same place
    // MeshView sometimes enable depth test
    glDisable(GL_DEPTH_TEST);

    /*
    if (redraw_all) {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

        // draw the image up-side down: y' = -y + height 
        // and fill the image to the viewport by setting proper scale.
        glTranslatef(0.0f, win_height, 0.0f);
        glPixelZoom(win_width * 1.0 / img_width, -win_height * 1.0 / img_height);

        glRasterPos2f(0.0f, 0.0f);    
        glDrawPixels(img_width, img_height, GL_RGB, GL_UNSIGNED_BYTE, byte_data);
    
        glFinish();
        redraw_all = false;
        return;
    }*/

    // TODO: no guarantee that the queue is not corrupted due to parallel
    //mtx.lock();
    //mtx.unlock();
    bool local_draw_block = num_blocks_in_queue > 0;
    if (local_draw_block) {
        //display_block->copy(block_queue_front());
        display_block = block_queue_front(); 
        if (! display_block) return;
    }
    
    if (local_draw_block) {
        glPixelZoom((Float)win_width / img_width, (Float)win_height / img_height);
        glTranslatef(display_block->x0, img_height - display_block->y0, 0.0f);
        
        glPixelZoom(1, -1);
        glRasterPos2f(0.f, 0.f);    // bottom left
        glDrawBuffer(GL_FRONT_AND_BACK);        
        glDrawPixels(display_block->x1 - display_block->x0, 
                     display_block->y1 - display_block->y0,
                     GL_RGB,
                     GL_UNSIGNED_BYTE,
                     display_block->data);

        /*
        // draw border
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_LINE_STRIP);
        glVertex2f(0.0f,                                    0.0f);
        glVertex2f(display_block->x1 - display_block->x0,   0.0f);
        glVertex2f(display_block->x1 - display_block->x0,   display_block->y0 - display_block->y1);
        glVertex2f(0.0f,                                    display_block->y0 - display_block->y1);
        glEnd();*/
    }

    if (local_draw_block) {
        block_queue_pop();
    }
}

void ImageBlockView::keyboard(unsigned char key, int x, int y) {
	switch(key) {
        case 's': 
        {
            char file[256];
            sprintf(file, tmp_file_template.c_str(), cur_frame);
            this->save_linear_last_complete_frame(file);
            break;
        }
        case 'x':
        {
            char file[256];
            sprintf(file, tmp_screen_template.c_str(), cur_frame);
            this->save_linear_screen(file);
            break;
        }
        case 't':
        {
            char file[256];
            sprintf(file, tmp_screen_template.c_str(), cur_frame);
            //this->save_tonemap_last_complete_frame(file);
            break;
        }
    }
}

void ImageBlockView::mouse(MouseButton button, MouseState state, int x, int y) {
    FrameBuffer *frame = scene->get_frame_buffer();
    ImageFloat *cur = frame->get_current_buffer();
    Rgb color((*cur)(y, x, 0), (*cur)(y, x, 1), (*cur)(y, x, 2));    
    Log::info() << "Frame      : (row, col, color) = (" << y << ", " << x << ", " << color << ")" << endn;

    ImageFloat *tm = frame->get_tonemap_buffer();
    Rgb tonemapped((*tm)(y, x, 0), (*tm)(y, x, 1), (*tm)(y, x, 2));    
    Log::info() << "Tonemapped : (row, col, color) = (" << y << ", " << x << ", " << tonemapped << ")" << endn;

    this->mouse_event.notify_click(button, state, x, y);
}

void ImageBlockView::motion(int x, int y) {    
}

void ImageBlockView::on_close() {
}

void ImageBlockView::save_linear_last_complete_frame(const char *file) {
    /*
    if (renderer->is_gpu_direct_output()) {                
        renderer->fill_progressive_buffer(scene, frame);

        frame->store_last_progressive_buffer();
    }*/
        
    ImageFloat img(*frame->get_last_complete_frame_buffer());
    
    if (img.save(file))
        Log::info() << file << " saved." << endn;
}

void ImageBlockView::save_linear_screen(const char *file) {
    /*
    if (renderer->is_gpu_direct_output()) {                
        renderer->fill_progressive_buffer(scene, frame);
    }*/
        
    ImageFloat img(*frame->get_current_buffer());
    
    if (img.save(file))
        Log::info() << file << " saved." << endn;
}

} // end namespace