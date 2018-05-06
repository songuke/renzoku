#include "frame.h"

namespace Renzoku {

FrameBuffer::FrameBuffer(Size2 img_size) {
    int height = img_size.height;
    int width  = img_size.width;

    this->byte_buffer               = new ImageByte((Byte)0, height, width, 3);
    this->tonemap_buffer            = new ImageFloat(0.0f, height, width, 3);        
    this->buffer                    = new ImageFloat(0.0f, height, width, 3);
    this->frame_buffer              = new ImageFloat(0.0f, height, width, 3);
    this->scale = 255.0f;

    this->index = 0;        
    this->display_correction_factor = 1;

    tonemapper = new ReinhardToneMap(0.18, img_size, 3);
    gamma_mapper = new GammaToneMap();
}

FrameBuffer::~FrameBuffer() {
    delete byte_buffer;
    delete tonemap_buffer;
    delete buffer;
    delete frame_buffer;        
}

void FrameBuffer::reset() {
    this->byte_buffer->set(0);
    this->tonemap_buffer->set(0);   
    this->buffer->set(0);
    this->frame_buffer->set(0);
    this->scale = 255.0f;
    this->index = 0;
    this->display_correction_factor = 1;
    tonemapper->reset();
}

void FrameBuffer::on_update_image(int tid, ImageFloat *img, int x0, int y0, int x1, int y1) {
    for (int i = y0; i < y1; ++i) {
        for (int j = x0; j < x1; ++j) {
            (*buffer)(i, j, 0) = (*img)(i, j, 0);
            (*buffer)(i, j, 1) = (*img)(i, j, 1);
            (*buffer)(i, j, 2) = (*img)(i, j, 2);
        }
    }

    if (index > 0) {
        tonemapper->map(buffer, tonemap_buffer, x0, y0, x1, y1);
    } else {
        gamma_mapper->map(buffer, tonemap_buffer, x0, y0, x1, y1);        
    }
        
    for (int i = y0; i < y1; ++i) {
        for (int j = x0; j < x1; ++j) {
            (*byte_buffer)(i, j, 0) = (Byte)((*tonemap_buffer)(i, j, 0) * scale);
            (*byte_buffer)(i, j, 1) = (Byte)((*tonemap_buffer)(i, j, 1) * scale);
            (*byte_buffer)(i, j, 2) = (Byte)((*tonemap_buffer)(i, j, 2) * scale);
        }
    }

    observable.notify_update_image(tid, byte_buffer, x0, y0, x1, y1);
}

void FrameBuffer::on_update_border(int tid, int x0, int y0, int x1, int y1) {
    // just forward to block view for highlight the border of the in-progress block
    observable.notify_update_border(tid, x0, y0, x1, y1);
}

void FrameBuffer::on_complete_frame(int tid, ImageFloat *img) {
    ++index;
    frame_buffer->copy(*buffer);
    
    tonemapper->precompute(frame_buffer);                       // for more accurate tone mapping

    observable.notify_complete_frame(tid, img);
        
    //scale = 255.0f / tonemap_buffer->max();
    scale = 255.0f;
    // Log::info() << "Scale: " << scale << endn;
}

} // end namespace