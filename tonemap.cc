#include "tonemap.h"
#include "image.h"

namespace Renzoku {

GammaToneMap::GammaToneMap(Float gamma) : gamma(gamma), inv_gamma(1.0f / gamma) {

}

void GammaToneMap::map(ImageFloat *src, ImageFloat *dst) {
    Size2 size = src->get_size();
    int height = size.height;
    int width  = size.width;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            (*dst)(i, j, 0) = clamp01(pow((*src)(i, j, 0), inv_gamma));
            (*dst)(i, j, 1) = clamp01(pow((*src)(i, j, 1), inv_gamma));
            (*dst)(i, j, 2) = clamp01(pow((*src)(i, j, 2), inv_gamma));
        }
    }
}

void GammaToneMap::map(ImageFloat *src, ImageFloat *dst, int x0, int y0, int x1, int y1) {
    for (int i = y0; i < y1; ++i) {
        for (int j = x0; j < x1; ++j) {
            (*dst)(i, j, 0) = clamp01(pow((*src)(i, j, 0), inv_gamma));
            (*dst)(i, j, 1) = clamp01(pow((*src)(i, j, 1), inv_gamma));
            (*dst)(i, j, 2) = clamp01(pow((*src)(i, j, 2), inv_gamma));
        }
    }
}

ReinhardToneMap::ReinhardToneMap(Size2 size, int color) : key(0.18f), buffer(NULL), scale(1.0f) {
    allocate_buffer(size, color);
}

ReinhardToneMap::ReinhardToneMap(Float key, Size2 size, int color) : key(key), buffer(NULL), scale(1.0f) {
    allocate_buffer(size, color);
}

void ReinhardToneMap::allocate_buffer(Size2 size, int color) {
    if (buffer) {
        if ((buffer->get_size() == size) && (buffer->get_color() == color))
            return;
        else 
            delete buffer;
    }
    buffer = new ImageFloat(0.0f, size.height, size.width, color);
}

void ReinhardToneMap::map(ImageFloat *src, ImageFloat *dst) {
    // Implement the simple equation 3 in Reinhard's SIGGRAPH '02 paper.
    if ((src->get_color() != 3) && (dst->get_color() != 3)) {
        Log::warn() << "Failed to tone map. Only three-color image is supported." << endn;
        return;
    }

    Size2 size = src->get_size();
    int height = size.height;
    int width  = size.width;

    Float *src_data = src->get_data();
    Float *dst_data = dst->get_data();
    Float *luminance = buffer->get_data();

    for (int i = 0; i < height * width; ++i) {
        if (is_nan(src_data[3 * i]) || is_nan(src_data[3 * i + 1]) || is_nan(src_data[3 * i + 2]))
            luminance[i] = 0;
        else
            luminance[i] = 0.27 * src_data[3 * i] + 0.67 * src_data[3 * i + 1] + 0.06 * src_data[3 * i + 2];                

        if (src_data[3 * i] < 0 || src_data[3 * i + 1] < 0 || src_data[3 * i + 2] < 0) {
            Log::info() << "Negative luminance. Tonemapped values not reliable." << endn;
        }
    }
    
    // find log-average luminance
    Float delta = 1e-4; // NOTE: change this parameter to control the saturation of the tone mapped image.
    Float avg_luminance = 0.;    
    Float prev_avg = 0;
    for (int i = 0; i < height * width; ++i) {        
        Float log = std::log(delta + luminance[i]);
        prev_avg = avg_luminance;
        avg_luminance += log;        
    }
    
    //cout << "sum log luminance: " << avg_luminance << endl;    
    avg_luminance = std::exp(avg_luminance / (height * width)); // equation 1
    //cout << "avg luminance: " << avg_luminance << endl;

    // scale 
    Float scale = key / avg_luminance;
    //cout << "scale: " << scale << endl;
    for (int i = 0; i < height * width; ++i) {
        Float org_luminance = luminance[i];
        luminance[i] *= scale;                              // equation 2
        luminance[i] = luminance[i] / (1. + luminance[i]);  // equation 3
        // trick: use luminance to store the ratio to scale RGB image
        if (org_luminance > 0)
            luminance[i] /= org_luminance;
        else
            luminance[i] = 0;
    }
    
    // new color = chrominance * new luminance 
    //           = old color / old luminance * new luminance
    //           = old color * the scale stored in the luminance array
    for (int i = 0; i < height * width; ++i) {
        dst_data[3 * i    ] = src_data[3 * i]     * luminance[i];
        dst_data[3 * i + 1] = src_data[3 * i + 1] * luminance[i];
        dst_data[3 * i + 2] = src_data[3 * i + 2] * luminance[i];

        dst_data[3 * i    ] = clamp01(dst_data[3 * i]);
        dst_data[3 * i + 1] = clamp01(dst_data[3 * i + 1]);
        dst_data[3 * i + 2] = clamp01(dst_data[3 * i + 2]);
    }
}

void ReinhardToneMap::precompute(ImageFloat *src) {
    // cache average luminance of the whole image
    if ((src->get_color() != 3)) {
        Log::warn() << "Failed to tone map. Only three-color image is supported." << endn;
        return;
    }

    Size2 size = src->get_size();
    int height = size.height;
    int width  = size.width;

    Float *src_data = src->get_data();
    Float *luminance = buffer->get_data();

    for (int i = 0; i < height * width; ++i) {
        if (is_nan(src_data[3 * i]) || is_nan(src_data[3 * i + 1]) || is_nan(src_data[3 * i + 2]))
            luminance[i] = 0;
        else
            luminance[i] = 0.27 * src_data[3 * i] + 0.67 * src_data[3 * i + 1] + 0.06 * src_data[3 * i + 2];                

        if (src_data[3 * i] < 0 || src_data[3 * i + 1] < 0 || src_data[3 * i + 2] < 0) {
            Log::info() << "Negative luminance. Tonemapped values not reliable." << endn;
        }
    }
    
    // find log-average luminance
    Float delta = 1e-4; // NOTE: change this parameter to control the saturation of the tone mapped image.
    Float avg_luminance = 0.;    
    Float prev_avg = 0;
    for (int i = 0; i < height * width; ++i) {        
        Float log = std::log(delta + luminance[i]);        
        prev_avg = avg_luminance;
        avg_luminance += log;        
    }
    
    //cout << "sum log luminance: " << avg_luminance << endl;    
    avg_luminance = std::exp(avg_luminance / (height * width)); // equation 1
    //cout << "avg luminance: " << avg_luminance << endl;

    scale = key / avg_luminance;    
}


void ReinhardToneMap::map(ImageFloat *src, ImageFloat *dst, int x0, int y0, int x1, int y1) {
    if ((src->get_color() != 3) && (dst->get_color() != 3)) {
        Log::warn() << "Failed to tone map. Only three-color image is supported." << endn;
        return;
    }

    int height = src->get_height();
    int width  = src->get_width();

    Float *src_data = src->get_data();
    Float *dst_data = dst->get_data();
    Float *luminance = buffer->get_data();

    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            int i = y * width + x;
    
            if (is_nan(src_data[3 * i]) || is_nan(src_data[3 * i + 1]) || is_nan(src_data[3 * i + 2]))
                luminance[i] = 0;   
            else
                luminance[i] = 0.27 * src_data[3 * i] + 0.67 * src_data[3 * i + 1] + 0.06 * src_data[3 * i + 2];                
            
            Float org_luminance = luminance[i];
            luminance[i] *= scale;                              // equation 2
            luminance[i] = luminance[i] / (1. + luminance[i]);  // equation 3
            
            // trick: use luminance to store the ratio to scale RGB image
            if (org_luminance > 0)
                luminance[i] /= org_luminance;
            else
                luminance[i] = 0;

            // new color = chrominance * new luminance 
            //           = old color / old luminance * new luminance
            //           = old color * the scale stored in the luminance array
            dst_data[3 * i    ] = src_data[3 * i]     * luminance[i];
            dst_data[3 * i + 1] = src_data[3 * i + 1] * luminance[i];
            dst_data[3 * i + 2] = src_data[3 * i + 2] * luminance[i];

            dst_data[3 * i    ] = clamp01(dst_data[3 * i]);
            dst_data[3 * i + 1] = clamp01(dst_data[3 * i + 1]);
            dst_data[3 * i + 2] = clamp01(dst_data[3 * i + 2]);
        }
    }
}

void ReinhardToneMap::reset() {
    buffer->set(0.f);
    scale = 1.0f;
}

} // end namespace