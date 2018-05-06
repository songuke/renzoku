#ifndef _SENSOR_H_
#define _SENSOR_H_

#include "math3.h"
#include "ray.h"

namespace Renzoku {

struct SensorSample {
	int x, y;
    Vec2 jitter;

    SensorSample() : jitter(Vec2(0.5f, 0.5f)) {
    }
};
typedef vector<SensorSample> SensorSamples;


struct SensorSamplePack {
    enum {
        SAMPLE_PACK_DEFAULT,
        SAMPLE_PACK_BEGIN_FRAME,
        SAMPLE_PACK_END_FRAME
    } tag;
	SensorSamples samples;
    int x0, y0, x1, y1;

    SensorSamplePack() : tag(SAMPLE_PACK_DEFAULT) {
    }
};
typedef vector<SensorSamplePack> SensorSamplePacks;

class Sensor;

class PixelAllocator {
public:
    PixelAllocator(Sensor *sensor);

    SensorSample    get_sensor_sample(Random &rd, Vec2 pixel);
    void            get_sensor_sample_2x2(Random &rd, Vec2 pixel, SensorSamples &samples);

private:
    Sensor *sensor;
};


class BlockAllocator {
public:
    BlockAllocator(PixelAllocator *pix, Sensor *sensor);
    ~BlockAllocator();

    void get_block(Random &rd, SensorSamplePack &pack);
    bool has_block() const;
    void reset();

private:
    PixelAllocator *pix;
    Sensor *sensor;

    int offset_x;
    int offset_y;
};


/**
 * A sensor contains a rectangular pixel collection. 
 *
 * The organization of samples in a pixel is specified by a PixelAllocator. 
 * How pixels are grouped into blocks are specified by BlockAllocator. 
 *
 * Sensor also handles anti-aliasing. It can be implemented by specifiying pixel position in the SensorRay class.
 * 
 * Different types of sensor: RGB, depth, norma, etc.
 * 
 */
class Sensor {
public:	
    Sensor();
    Sensor(Size2 img_size, Size2 film_size);
    
    void reset();
    void get_sample_packs(int num_packs, Random &rd, SensorSamplePacks &packs, bool &has_end_frame);
    
    inline Size2 get_block_size() const;
    inline void  set_block_size(Size2 block_size);
    
    inline Size2 get_image_size() const;
    inline void  set_image_size(Size2 img_size);

    inline Size2 get_film_size() const;
    inline void  set_film_size(Size2 film_size);
    
    inline Vec2 get_cop() const;
    inline Vec2 get_image_film_ratio() const;

    /**
     * Only render a set of small pixels for debugging purpose/collecting data.
     */
    inline void set_debug_pixels(const vector<Vec2> &pixels);
    inline bool has_debug_pixels() const;
    void get_debug_sample_packs(SensorSamplePacks &packs);

private:
	PixelAllocator *pixel_alloc;	
	BlockAllocator *block_alloc;
    
    Size2 block_size;    
    Size2 film_size; 
    Size2 img_size;

    vector<Vec2> debug_pixels;
};


enum SensorType {
    SENSOR_RGB,
    SENSOR_DIRECT,
    SENSOR_INDIRECT,
    SENSOR_DEPTH,
    SENSOR_NORMAL
};


struct PrimaryRay {
    Ray r;
    SensorSample s;
    SensorType type;
    Rgb val;

    PrimaryRay(const Ray &r, const SensorSample &s) 
        : r(r), s(s), type(SENSOR_RGB), val(DefaultRgb::black) {
    }
};
typedef vector<PrimaryRay> PrimaryRays;


struct PrimaryRayBundle {
    enum RayBundleEnum {
        BUNDLE_DEFAULT,
        BUNDLE_BEGIN_FRAME,
        BUNDLE_END_FRAME
    } tag;

    PrimaryRays rays;
    int x0, y0, x1, y1;     // image block to update

    PrimaryRayBundle() : tag(BUNDLE_DEFAULT) {
    }
};
typedef vector<PrimaryRayBundle> PrimaryRayBundles;

inline Size2 Sensor::get_block_size() const {
    return block_size;
}

inline Size2 Sensor::get_film_size() const {
    return film_size;
}

inline Size2 Sensor::get_image_size() const {
    return img_size;
}

inline void Sensor::set_block_size(Size2 block_size) {
    this->block_size = block_size;
}

inline void Sensor::set_image_size(Size2 image_size) {
    this->img_size = image_size;
}

inline void Sensor::set_film_size(Size2 film_size) {
    this->film_size = film_size;
}

inline Vec2 Sensor::get_cop() const {        
    return (img_size / 2.0f).to_vec2();
}

inline Vec2 Sensor::get_image_film_ratio() const {
    return (img_size / film_size).to_vec2();
}

inline void Sensor::set_debug_pixels(const vector<Vec2> &pixels) {
    debug_pixels = vector<Vec2>(pixels.begin(), pixels.end());
}

inline bool Sensor::has_debug_pixels() const {
    return debug_pixels.size() > 0;
}

} // end namespace

#endif