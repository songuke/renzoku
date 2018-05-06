#include "sensor.h"
#include "log.h"

namespace Renzoku {

PixelAllocator::PixelAllocator(Sensor *sensor) : sensor(sensor) {
}

SensorSample PixelAllocator::get_sensor_sample(Random &rd, Vec2 pixel) {
    SensorSample s;
    // Inaccuracy can occur: 511 + rd() 
    //   where rd() returns 0x00ffffff / (1 + 0x00ffffff) = 0.999999940395355224609375
    //   = 0.99999994f = 0x3f7fffff in 32-bit float.
    // The sum 511.99999994 will be rounded to 512 in 32-bit float.
    // This will cause index out of bound when color is written to the pixel in the output image (pixel 512 not exist).

    // I choose to store the integer coordinates of the pixel separately.
    // This is redundant but ensure no out of bound error occurs.
    s.x = pixel.x();
    s.y = pixel.y();
    s.jitter = Vec2(rd(), rd());
    return s;
}

void PixelAllocator::get_sensor_sample_2x2(Random &rd, Vec2 pixel, SensorSamples &samples) {
    // 2x2 stratified sampling
    SensorSample s;
    s.x = pixel.x();
    s.y = pixel.y();
    s.jitter = Vec2(rd() * 0.5f,        rd() * 0.5f);
    samples.push_back(s);

    s.jitter = Vec2(rd() * 0.5f + 0.5f, rd() * 0.5f);
    samples.push_back(s);
    
    s.jitter = Vec2(rd() * 0.5f,        rd() * 0.5f + 0.5f);
    samples.push_back(s);
    
    s.jitter = Vec2(rd() * 0.5f + 0.5f, rd() * 0.5f + 0.5f);
    samples.push_back(s);
}

BlockAllocator::BlockAllocator(PixelAllocator *pix, Sensor *sensor) 
    : pix(pix), sensor(sensor)
{    
    reset();
}

BlockAllocator::~BlockAllocator() { 
}

void BlockAllocator::get_block(Random &rd, SensorSamplePack &pack) {
    Size2 block_size = sensor->get_block_size();
    Size2 image_size = sensor->get_image_size();

    for (int i = 0; i < block_size.height; ++i) {
        for (int j = 0; j < block_size.width; ++j) {
            Vec2 pixel(offset_x + j, offset_y + i);
            //pack.samples.push_back(pix->get_sensor_sample(rd, pixel));
            pix->get_sensor_sample_2x2(rd, pixel, pack.samples);
        }
    }
    pack.x0 = offset_x;
    pack.x1 = offset_x + block_size.width;
    pack.y0 = offset_y;
    pack.y1 = offset_y + block_size.height;
        
    offset_x += block_size.width;   
    if (offset_x >= image_size.width) {
        offset_x = 0;
        offset_y += block_size.height;        
    }

    if (offset_y >= image_size.height) {        
        pack.tag = SensorSamplePack::SAMPLE_PACK_END_FRAME;
        offset_x = 0;
        offset_y = 0;
    }
}

void BlockAllocator::reset() {
    offset_x = 0;
    offset_y = 0;
}

Sensor::Sensor() 
    : img_size(Size2(512, 512)), film_size(0.025, 0.025) {
    block_size = Size2(RAY_BUNDLE_DIM, RAY_BUNDLE_DIM);
    pixel_alloc = new PixelAllocator(this);
    block_alloc = new BlockAllocator(pixel_alloc, this);
}

Sensor::Sensor(Size2 img_size, Size2 film_size) 
    : img_size(img_size), film_size(film_size) {
    block_size = Size2(RAY_BUNDLE_DIM, RAY_BUNDLE_DIM);
    pixel_alloc = new PixelAllocator(this);
    block_alloc = new BlockAllocator(pixel_alloc, this);
}

void Sensor::reset() {
    block_alloc->reset();
}

void Sensor::get_sample_packs(int num_packs, Random &rd, SensorSamplePacks &packs, bool &has_end_frame) {
    has_end_frame = false;
    while (num_packs > 0) {
        SensorSamplePack pack;
        block_alloc->get_block(rd, pack);
        packs.push_back(pack);

        if (pack.tag == SensorSamplePack::SAMPLE_PACK_END_FRAME)
            has_end_frame = true;

        --num_packs;
    }
}

void Sensor::get_debug_sample_packs(SensorSamplePacks &packs) {    
    for (int i = 0; i < debug_pixels.size(); ++i) {        
        SensorSample s;    
        s.x = debug_pixels[i].x();
        s.y = debug_pixels[i].y();
        s.jitter = Vec2(0, 0);

        SensorSamplePack pack;
        pack.samples.push_back(s);
        pack.x0 = s.x;
        pack.y0 = s.y;
        pack.x1 = s.x + 1;
        pack.y1 = s.y + 1;
        packs.push_back(pack);
    }
}

} // end namespace