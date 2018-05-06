#ifndef _TONEMAP_H_
#define _TONEMAP_H_

#include "common.h"
#include "math3.h"

namespace Renzoku {

class ToneMap {
public:
    virtual void map(ImageFloat *src, ImageFloat *dst) = 0;
    virtual void map(ImageFloat *src, ImageFloat *dst, int x0, int y0, int x1, int y1) = 0;

    inline static Float clamp01(Float a) {
        return a > 1.f ? 1.f : (a < 0.f ? 0.f : a);
    }

    virtual void reset() {}
};

class GammaToneMap : public ToneMap {
public:
    GammaToneMap(Float gamma = 2.2f);

    virtual void map(ImageFloat *src, ImageFloat *dst);
    virtual void map(ImageFloat *src, ImageFloat *dst, int x0, int y0, int x1, int y1);

protected:
    Float gamma;
    Float inv_gamma;
};

class ReinhardToneMap : public ToneMap {
public:
    ReinhardToneMap(Size2 size, int color);
    ReinhardToneMap(Float key, Size2 size, int color);

    /**
     * Tone map the image from src to dst.
     */
    virtual void map(ImageFloat *src, ImageFloat *dst);

    /**
     * Calculate the scale value from a source image. 
     * 
     * The scale value is used for tone mapping of blocks.
     */
    void precompute(ImageFloat *src);

    /**
     * Tone map a region with precomputed parameters.
     * 
     * precompute_tonemap(src) must be called before.
     */ 
    virtual void map(ImageFloat *src, ImageFloat *dst, int x0, int y0, int x1, int y1);

    virtual void reset();

private:
    void allocate_buffer(Size2 size, int color);

private:
    Float key;
    Float scale;            // parameter for tone mapping
    ImageFloat *buffer;
};

} // end namespace

#endif