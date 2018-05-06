#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <cstdio>
#include <cstring>
#include <typeinfo>
#include <cmath>
#include <cfloat>
using namespace std;

//#ifdef _WIN32
//#define hypotf hypot // For the OpenEXR headers
//#endif
#include <ImfInputFile.h>
#include <ImfRgbaFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>

#include "common.h"

#ifdef FREEIMAGE
#ifdef _MSC_VER
#include <Windows.h>    // add Windows.h before FreeImage to avoid BOOL to be defined by FreeImage
#endif
#include <FreeImage.h>
#endif

#include "error.h"
#include "vec2.h"
#include "size.h"
#include "rgb.h"
#include "log.h"

namespace Renzoku {
template <class T>
class Image {
 public:
	Image();
	Image(const Image &);
	Image(T val, int height, int width);
	Image(T val, int height, int width, int color);
	Image(T *buf, int height, int width);
	Image(T *buf, int height, int width, int color);
	Image(const char *file);
	~Image();

    void initialize(T val, int height, int width);
	void initialize(T val, int height, int width, int color);	

 public:
    int get_height() const;
    int get_width() const;
    int get_color() const;
    Size2 get_size() const;

 public:
	T operator[](int i) const; // 1D access
	T& operator[](int i);
	T operator()(int i, int j) const; // 2D gray access
	T& operator()(int i, int j);
	T operator()(int i, int j, int k) const; // 2D color access
	T& operator()(int i, int j, int k);

    Rgb lookup(Float u, Float v) const;

 public:
    bool save(const char *file);
	
	void load(const char *file);
	void load(T *buf, int height, int width);
	void load(T *buf, int height, int width, int color);

	void copy(const Image<T> &img);

    template <class T2>
    void copy(const Image<T2> &img);

 public:
	T *get_data() const;
 public:
	void to_gray();		

	void linearize(T *inv_response);
	void linearize(T *inv_r, T *inv_g, T *inv_b, float log_exposure_time);	
	// linearize and normalize to [0, 255]
	void linearize_8bit(T *inv_r, T *inv_g, T *inv_b, float log_exposure_time);

	void add(const Image<T> &img);
	void subtract(const Image<T> &img);
    void multiply(const Image<T> &img);
    void divide(const Image<T> &img);    
    void diff(const Image<T> &img, const Image<T> &diff) const;

	void add(T val);
	void multiply(T val);

    void set(T val);
	void set(const Vec2 &pixel, const Rgb &val);
	void accumulate(const Vec2 &pixel, const Rgb &val);
	bool is_valid_pixel(const Vec2 &pixel) const;

	void log();
	void exp();
    void tonemap(Float key = 0.18);
    void tonemap(int x0, int y0, int x1, int y1, Float key = 0.18);
	

    /**
     * Clamp pixel values to [0, 1]
     */
    void clamp(T new_min = 0, T new_max = 1);

	T max();
	T min();

	void normalize();
	void normalize_to(T new_min = 0, T new_max = 1);

    void draw_line(const Vec2 &from, const Vec2 &to, const Rgb &color, Float alpha);

    void flip_vertical();

protected:
    bool save_ppm(const char *file);
	bool save_exr(const char *file);

private:
	T *data;
	
    int height, width;
	int color;
	int bit_depth;
};


template <class T>
Image<T>::Image() : data(NULL), height(0), width(0), color(0), bit_depth(0) {
}

template <class T>
Image<T>::Image(const Image &img) {
	height = img.height;
	width = img.width;
	color = img.color;
	bit_depth = 8;
	if (img.data) {
		data = new T[height * width * color];
		if (!data) error_alloc(__FILE__, __LINE__);
		memcpy(data, img.data, sizeof(T) * height * width * color);
	}
}

template <class T>
Image<T>::Image(T val, int height, int width) {
	data = new T[height * width];
	if (!data) error_alloc(__FILE__, __LINE__);
	this->height = height;
	this->width = width;	
	this->color = 1;
	this->bit_depth = 8;
	for (int i = 0; i < height * width; ++i)
		data[i] = val;
}

template <class T>
Image<T>::Image(T val, int height, int width, int color) {
	data = new T[height * width * color];
	if (!data) error_alloc(__FILE__, __LINE__);
	this->height = height;
	this->width = width;	
	this->color = color;
	this->bit_depth = 8;
	for (int i = 0; i < height * width * color; ++i)
		data[i] = val;
}

template <class T>
Image<T>::Image(T *buf, int height, int width) {
	data = new T[height * width];
	if (!data) error_alloc(__FILE__, __LINE__);
	this->height = height;
	this->width = width;	
	this->color = 1;
	this->bit_depth = 8;
	memcpy(data, buf, sizeof(T) * height * width);
}

template <class T>
Image<T>::Image(T *buf, int height, int width, int color) {
	data = new T[height * width * color];
	if (!data) error_alloc(__FILE__, __LINE__);
	this->height = height;
	this->width = width;	
	this->color = color;
	this->bit_depth = 8;
	memcpy(data, buf, sizeof(T) * height * width * color);
}

template <class T>
Image<T>::Image(const char *file) {
	load(file);
}

template <class T>
Image<T>::~Image() {
	if (data)
		delete [] data;
}

template <class T>
void Image<T>::initialize(T val, int height, int width) {
    if (data && (height * width > this->height * this->width)) {
        delete [] data;
        data = NULL;
    }

    this->height = height;
	this->width = width;	
	this->color = 1;
	this->bit_depth = 8;

    if (!data)
        data = new T[this->height * this->width * this->color];
    if (!data) error_alloc(__FILE__, __LINE__);

	for (int i = 0; i < height * width * color; ++i)
		data[i] = val;
}

template <class T>
void Image<T>::initialize(T val, int height, int width, int color) {
    if (data && (height * width * color > this->height * this->width * this->color)) {
        delete[] data;
        data = NULL;
    }
    
	this->height = height;
	this->width = width;	
	this->color = color;
	this->bit_depth = 8;
	
    if (!data)
        data = new T[this->height * this->width * this->color];
    if (!data) error_alloc(__FILE__, __LINE__);

    for (int i = 0; i < height * width * color; ++i)
		data[i] = val;
}


template <class T>
int Image<T>::get_height() const { 
    return height; 
}

template <class T>
int Image<T>::get_width() const { 
    return width; 
}

template <class T>
int Image<T>::get_color() const {
    return color;
}

template <class T>
Size2 Image<T>::get_size() const {
    return Size2(width, height);
}

template <class T>
T Image<T>::operator[](int i) const {
	return data[i];
}

template <class T>
T &Image<T>::operator[](int i) {
	return data[i];
}

template <class T>
T Image<T>::operator()(int i, int j) const {
	return data[i * width + j];
}

template <class T>
T &Image<T>::operator()(int i, int j) {
	return data[i * width + j];
}

template <class T>
T Image<T>::operator()(int i, int j, int k) const {
	return data[color * (i * width + j) + k];
}

template <class T>
T &Image<T>::operator()(int i, int j, int k) {
	return data[color * (i * width + j) + k];
}

template <class T>
bool Image<T>::save_ppm(const char *file) {
	FILE *f = fopen(file, "wb");
	if (!f) error_file(__FILE__, __LINE__, file);
	
	switch (color) {
		case 1:
			if (typeid(T) == typeid(Byte)) {
				fprintf(f, "P5\n%d %d\n255\n", width, height);
				fwrite(data, sizeof(Byte), height * width, f);
			} else {
				// TODO: should save as PFM
				throw_warning("saving image as 8-bit format; accuracy lost.", __FILE__, __LINE__);
				fprintf(f, "P5\n%d %d\n255\n", width, height);
				Byte *buf = new Byte[height * width];
				for (int i = 0; i < height * width; ++i)
					buf[i] = (Byte)data[i];
				fwrite(buf, sizeof(Byte), height * width, f);
				delete [] buf;
			}
			break;
		case 3:
			if (typeid(T) == typeid(Byte)) {
				fprintf(f, "P6\n%d %d\n255\n", width, height);
				fwrite(data, sizeof(Byte), height * width * color, f);
			} else {
				// TODO: should save as PFM
				throw_warning("saving image as 8-bit format; accuracy lost.", __FILE__, __LINE__);
				fprintf(f, "P6\n%d %d\n255\n", width, height);
				Byte *buf = new Byte[height * width * color];
				for (int i = 0; i < height * width * color; ++i)
					buf[i] = (Byte)data[i];
				fwrite(buf, sizeof(Byte), height * width * color, f);
				delete [] buf;
			}
			break;
	}
	fclose(f);
    return true;
}

template <class T>
bool Image<T>::save_exr(const char *filename) {
	if (typeid(T) != typeid(Float) && typeid(T) != typeid(int)) {
		throw_error("only float/int image type is currently supported.", 
					__FILE__, __LINE__);
	}

	if (color != 1 && color != 3) {
		throw_error("only grayscale or RGB color image is currently supported.", 
					__FILE__, __LINE__);
	}

    // NOTE: both overflow/NaN case will be displayed as 1.#J when the image is viewed with exrdisplay.
    // check for underflow and overflow
    bool warning = false;
    //const Float EXR_MIN = 6.1e-5f; // from half.h in OpenEXR
    const Float EXR_MAX = 6.5e+4f;
    for (int i = 0; i < height * width * color; ++i) {
        if (data[i] > EXR_MAX) {
            warning = true;
            break;
        }
    }
    if (warning) {
        Log::warn() << "Overflow during conversion from 32-bit (or 64-bit) floating point to OpenEXR half (16-bit)." << '\n';
        Log::warn() << "Image min/max: " << min() << "\t" << max() << '\n';
    }

	// call OpenEXR function to save
	Imf::Rgba *pixels = new Imf::Rgba[height * width];
	for (int i = 0; i < height * width; ++i) {
        if (color == 3) {
		    pixels[i] = Imf::Rgba(data[3 * i], data[3 * i + 1], data[3 * i + 2], 1.0f);
        } else if (color == 1) {
            pixels[i] = Imf::Rgba(data[i], data[i], data[i], 1.0f);
        }
	}

	Imf::RgbaOutputFile file (filename, width, height, Imf::RgbaChannels::WRITE_RGBA);
	file.setFrameBuffer (pixels, 1, width);
	file.writePixels (height);
	
	/*
	Header header (width, height);
	header.channels().insert ("Z", Channel (FLOAT));
	OutputFile file (filename, header); 
	FrameBuffer frameBuffer; 

	frameBuffer.insert ("Z", 
						Slice (FLOAT, 
							   (char *) data, 
							   sizeof (*data) * 1, 
							   sizeof (*data) * width)); 
	file.setFrameBuffer (frameBuffer);
	file.writePixels (height);
	*/

    delete[] pixels;

    return true;
}

template <class T>
bool Image<T>::save(const char *file) {
    // determine extension and forward to proper save functions
    char *ext = (char *)file + strlen(file);
    do {
        --ext;
    } while (ext >= file && *ext != '.');
	++ext;

    if (ext < file) {
        Log::warn() << "Failed to save image. Unknown image format." << endn;
        return false;
    }
    else if (strcmp(ext, "ppm") == 0) {
        return save_ppm(file);
    }    
    else if (strcmp(ext, "exr") == 0) {
        return save_exr(file);
    } else {
        Log::warn() << "Failed to save image. Unknown image format." << endn;
		return false;
	}
    return true;
}

template <class T>
void Image<T>::load(const char *file) {
	FILE *f = fopen(file, "rb");
	if (!f) error_file(__FILE__, __LINE__, file);
	
	char buf[8];
	int maxval;
	fgets(buf, 8, f); // get header
	if (buf[0] == 'P' && buf[1] == '5') { // PGM
		fscanf(f, "%d %d\n", &width, &height);
		fscanf(f, "%d\n", &maxval);
		
		bit_depth = 8; // TODO: support 16-bit PPM
		color = 1;
		data = new T[height * width];
		if (!data) error_alloc(__FILE__, __LINE__);

		if (typeid(T) == typeid(Byte)) {
			fread(data, sizeof(Byte), height * width, f);
		} else {
			Byte *buf = new Byte[height * width];
			if (!buf) error_alloc(__FILE__, __LINE__);
			fread(buf, sizeof(Byte), height * width, f);
				
			for (int i = 0; i < height * width; ++i)
				data[i] = (T)buf[i] / 255.0f;

			delete [] buf;
		}

        fclose(f);
	} else if (buf[0] == 'P' && buf[1] == '6') { // PPM
		fscanf(f, "%d %d\n", &width, &height);
		fscanf(f, "%d\n", &maxval);

		bit_depth = 8;
		color = 3;
		data = new T[height * width * color];
		if (!data) error_alloc(__FILE__, __LINE__);

		if (typeid(T) == typeid(Byte)) {
			fread(data, sizeof(Byte), height * width * color, f);
		} else {
			Byte *buf = new Byte[height * width * color];
			if (!buf) error_alloc(__FILE__, __LINE__);
			fread(buf, sizeof(Byte), height * width * color, f);
				
			for (int i = 0; i < height * width * color; ++i)
				data[i] = (T)buf[i] / 255.0f;

			delete [] buf;
		}

        fclose(f);

	} else {

        fclose(f);

#ifdef MAGICK
		// call Magick++ to handle other files
		Magick::Image img; // create an *empty* image using the default Image constructor
		img.read(file);
		
		height = img.rows();
        width = img.columns();
		bit_depth = img.modulusDepth();
		color = 3;
		data = new T[height * width * color];
		
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {				
                Magick::ColorRGB color(img.pixelColor(j, i));
				data[(i * width + j) * 3    ] = (T)color.red();
				data[(i * width + j) * 3 + 1] = (T)color.green();
				data[(i * width + j) * 3 + 2] = (T)color.blue();
            }
        }

#else
#ifdef FREEIMAGE

        // check the file signature and deduce its format
        // (the second argument is currently not used by FreeImage)
        FREE_IMAGE_FORMAT fif = FreeImage_GetFileType(file, 0);
        if (fif == FIF_UNKNOWN) {
            // no signature, try to guess the file format from the file extension
            fif = FreeImage_GetFIFFromFilename(file);
        }
        
        // check that the plugin has reading capabilities ...
        if ((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
            Log::info() << "Loading image : " << file << endn;

            // OK, let's load the file  
            FIBITMAP *dib = FreeImage_Load(fif, file, 0);
            
            // calculate the number of bytes per pixel (3 for 24-bit or 4 for 32-bit)

            FREE_IMAGE_TYPE image_type = FreeImage_GetImageType(dib);

            BYTE *bits = FreeImage_GetBits(dib);
            unsigned pitch = FreeImage_GetPitch(dib);
            width = FreeImage_GetWidth(dib);
            height = FreeImage_GetHeight(dib);

            int bytespp = FreeImage_GetLine(dib) / FreeImage_GetWidth(dib);
            bit_depth = bytespp * 8;		

            // always assume 3 channels
            color = 3;
            data = new T[height * width * color];

            if (image_type == FIT_BITMAP) {
                for(unsigned y = 0; y < FreeImage_GetHeight(dib); y++) {
                    BYTE *pixel = (BYTE *)bits; // beginning of a line
                    for(unsigned x = 0; x < FreeImage_GetWidth(dib); x++) {                 
                        data[3 * (y * width + x)    ] = pixel[FI_RGBA_RED] / 255.f;
                        data[3 * (y * width + x) + 1] = pixel[FI_RGBA_GREEN] / 255.f;
                        data[3 * (y * width + x) + 2] = pixel[FI_RGBA_BLUE] / 255.f;

                        // jump to next pixel
                        pixel += bytespp;
                    }
                    bits += pitch;
                }
            } else if (image_type == FIT_RGBF) {
                for(unsigned y = 0; y < FreeImage_GetHeight(dib); y++) {
                    FIRGBF *pixel = (FIRGBF *)bits; // beginning of a line
                    for(unsigned x = 0; x < FreeImage_GetWidth(dib); x++) {                 
                        data[3 * (y * width + x)    ] = pixel[x].red;
                        data[3 * (y * width + x) + 1] = pixel[x].green;
                        data[3 * (y * width + x) + 2] = pixel[x].blue;
                    }
                    bits += pitch;
                }
            } else {
                    Log::warn() << "FreeImage: unimplemented image type." << endn;
            }            
        } else {
            Log::warn() << "FreeImage: image format not supported: " << file << endn;       
        }

#else
		throw_error("not a PGM or PPM file.", __FILE__, __LINE__);
#endif
#endif
	}
}

template <class T>
void Image<T>::load(T *buf, int height, int width) {
	if (!data) {
		data = new T[height * width];
		if (!data) error_alloc(__FILE__, __LINE__);

		this->color = 1;
		this->height = height;
		this->width = width;
	}
	memcpy(data, buf, sizeof(T) * height * width);
}

template <class T>
void Image<T>::load(T *buf, int height, int width, int color) {
	if (!data) {
		data = new T[height * width * color];
		if (!data) error_alloc(__FILE__, __LINE__);

		this->color = color;
		this->height = height;
		this->width = width;
	}
	memcpy(data, buf, sizeof(T) * height * width * color);
}

template <class T>
void Image<T>::copy(const Image<T>& img) {
    if (img.data) {
        bool alloc = false;
        if (!data) {
            alloc = true;
        } else {
            if (height != img.height || width != img.width || color != img.color) {
                delete [] data;
                alloc = true;
            }
        }
        if (alloc) data = new T[img.height * img.width * img.color];
		if (!data) error_alloc(__FILE__, __LINE__);
		memcpy(data, img.data, sizeof(T) * img.height * img.width * img.color);
	}
	height = img.height;
	width = img.width;
	color = img.color;
}

template <class T1>
template <class T2>
void Image<T1>::copy(const Image<T2>& img) {
    T2* img_data = img.get_data();
    int img_height = img.get_height();
    int img_width  = img.get_width();
    int img_color  = img.get_color();

    if (img_data) {
        bool alloc = false;
        if (!data) {
            alloc = true;
        } else {
            if (height != img_height || width != img_width || color != img_color) {
                delete [] data;
                alloc = true;
            }
        }
        if (alloc) data = new T1[img_height * img_width * img_color];
		if (!data) error_alloc(__FILE__, __LINE__);
		
        for (int i = 0; i < img_height * img_width * img_color; ++i)
            data[i] = (T1)img_data[i];
	}
	height = img_height;
	width  = img_width;
	color  = img_color;
}

template <class T>
T *Image<T>::get_data() const {
	return data;
}

template <class T>
void Image<T>::to_gray() {
	if (color == 1) return;

	if (typeid(T) == typeid(Byte)) {
		unsigned int sum; // avoid overflow
		for (int i = 0; i < height * width; ++i) {
			sum = 0;
			for (int k = 0; k < color; ++k) {
				sum += data[color * i + k];
			}
			data[i] = sum / color;
		}
		color = 1;
	} else {
		T sum;
		for (int i = 0; i < height * width; ++i) {
			sum = 0;
			for (int k = 0; k < color; ++k) {
				sum += data[color * i + k];
			}
			data[i] = sum / color;
		}
		color = 1;
	}
}

template <class T>
void Image<T>::add(const Image<T> &img) {
	if (typeid(T) == typeid(Byte)) {
		int val;
		for (int i = 0; i < height * width * color; ++i) {
			val = (int)data[i] + (int)img.data[i];
			if (val < 0) val = 0;
			if (val > 255) val = 255;
			data[i] = (Byte)val;
		}
	} else {
		for (int i = 0; i < height * width * color; ++i) {			
			data[i] = data[i] + img.data[i];
		}
	}
}

template <class T>
void Image<T>::subtract(const Image<T> &img) {
	if (typeid(T) == typeid(Byte)) {
		int val;
		for (int i = 0; i < height * width * color; ++i) {
			val = (int)data[i] - (int)img.data[i];
			if (val < 0) val = 0;
			if (val > 255) val = 255;
			data[i] = (Byte)val;
		}
	} else {		
		for (int i = 0; i < height * width * color; ++i) {			
			data[i] = data[i] - img.data[i];
		}
	}
}

template <class T>
void Image<T>::diff(const Image<T> &img, const Image<T> &diff) const {
	if (typeid(T) == typeid(Byte)) {
		int val;
		for (int i = 0; i < height * width * color; ++i) {
			val = (int)data[i] - (int)img.data[i];
			if (val < 0) val = 0;
			if (val > 255) val = 255;
			diff.data[i] = (Byte)val;
		}
	} else {		
		for (int i = 0; i < height * width * color; ++i) {			
			diff.data[i] = this->data[i] - img.data[i];
		}
	}
}

template <class T>
void Image<T>::multiply(const Image<T> &img) {
	if (typeid(T) == typeid(Byte)) {
		int val;
		for (int i = 0; i < height * width * color; ++i) {
			val = (int)data[i] * (int)img.data[i];
			if (val < 0) val = 0;
			if (val > 255) val = 255;
			data[i] = (Byte)val;
		}
	} else {		
		for (int i = 0; i < height * width * color; ++i) {			
			data[i] = data[i] * img.data[i];
		}
	}
}

template <class T>
void Image<T>::divide(const Image<T> &img) {
	if (typeid(T) == typeid(Byte)) {
		int val;
		for (int i = 0; i < height * width * color; ++i) {
            if (img.data[i] != 0) { // avoid division by zero
			    val = (int)data[i] / (int)img.data[i];
			    if (val < 0) val = 0;
			    if (val > 255) val = 255;            
            }
			data[i] = (Byte)val;
		}
	} else {		
		for (int i = 0; i < height * width * color; ++i) {			
            if (img.data[i] != 0) { // avoid division by zero
			    data[i] = data[i] / img.data[i];
            }
		}
	}
}

template <class T>
T Image<T>::max() {
	T max = FLT_MIN;
	for (int i = 0; i < height * width * color; ++i) 
		if (data[i] > max) max = data[i];
	return max;
}

template <class T>
T Image<T>::min() {
	T min = FLT_MAX;
	for (int i = 0; i < height * width * color; ++i) 
		if (data[i] < min) min = data[i];
	return min;
}

template <class T>
void Image<T>::linearize(T *inv_response) {
	if (color == 1) {
		for (int i = 0; i < height * width; ++i) {
			data[i] = inv_response[data[i]];
		}
	} else {
		throw_error("Invalid response function for color image.\n", __FILE__, __LINE__);
	}
}

template <class T>
void Image<T>::linearize(T *inv_r, T *inv_g, T *inv_b, float log_exposure_time) {	
	if (color == 3) {
		for (int i = 0; i < height * width; ++i) {
			float r = inv_r[(int)data[3 * i    ]];
			float g = inv_g[(int)data[3 * i + 1]];
			float b = inv_b[(int)data[3 * i + 2]];
			data[3 * i    ] = std::exp(r - log_exposure_time);
			data[3 * i + 1] = std::exp(g - log_exposure_time);
			data[3 * i + 2] = std::exp(b - log_exposure_time);
		}
	} else {
		throw_error("Invalid response function for gray image.\n", __FILE__, __LINE__);
	}
}

template <class T>
void Image<T>::linearize_8bit(T *inv_r, T *inv_g, T *inv_b, float log_exposure_time) {		if (color == 3) {
		/* map min_resp to 0 and max_resp to 255 */
		float min_r = FLT_MAX;
		float min_g = FLT_MAX;
		float min_b = FLT_MAX;
		float max_r = FLT_MIN;
		float max_g = FLT_MIN;
		float max_b = FLT_MIN; 
		for (int i = 0; i < 256; ++i) {
			if (inv_r[i] < min_r) min_r = inv_r[i];
			if (inv_r[i] > max_r) max_r = inv_r[i];
			if (inv_g[i] < min_g) min_g = inv_g[i];
			if (inv_g[i] > max_g) max_g = inv_g[i];
			if (inv_b[i] < min_b) min_b = inv_b[i];
			if (inv_b[i] > max_b) max_b = inv_b[i];
		}
		min_r = std::exp(min_r - log_exposure_time);
		max_r = std::exp(max_r - log_exposure_time);
		min_g = std::exp(min_g - log_exposure_time);
		max_g = std::exp(max_g - log_exposure_time);
		min_b = std::exp(min_b - log_exposure_time);
		max_b = std::exp(max_b - log_exposure_time);
				
		for (int i = 0; i < height * width; ++i) {
			float r = inv_r[(int)data[3 * i    ]];
			float g = inv_g[(int)data[3 * i + 1]];
			float b = inv_b[(int)data[3 * i + 2]];
			data[3 * i    ] = std::exp(r - log_exposure_time);
			data[3 * i + 1] = std::exp(g - log_exposure_time);
			data[3 * i + 2] = std::exp(b - log_exposure_time);
			data[3 * i    ] = (data[3 * i    ] - min_r) * 255 / (max_r - min_r);
			data[3 * i + 1] = (data[3 * i + 1] - min_g) * 255 / (max_g - min_g);
			data[3 * i + 2] = (data[3 * i + 2] - min_b) * 255 / (max_b - min_b);
		}
	} else {
		throw_error("Invalid response function for gray image.\n", __FILE__, __LINE__);
	}
}

template <class T>
void Image<T>::normalize() {
	if (typeid(T) != typeid(Byte)) { 
		T max_val = (1 << bit_depth) - 1;
		for (int i = 0; i < height * width * color; ++i)
			data[i] /= max_val;
	} else {
		throw_error("method not implemented.", __FILE__, __LINE__);
	}
}

template <class T>
void Image<T>::normalize_to(T new_min, T new_max) {
	/* scale current min and max to new_min and new_max */
	T max_val = max();
	T min_val = min();
	if (typeid(T) != typeid(Byte)) { 
		for (int i = 0; i < height * width * color; ++i)
			data[i] = ((data[i] - min_val) / (max_val - min_val)) * 
				(new_max - new_min) + new_min;
	} else {
		throw_error("method not implemented.", __FILE__, __LINE__);
	}
}


template <class T>
void Image<T>::clamp(T new_min, T new_max) {
    for (int i = 0; i < height * width * color; ++i) {
		if (data[i] < new_min) data[i] = new_min;
        if (data[i] > new_max) data[i] = new_max;
    }
}

template <class T>
void Image<T>::add(T val) {
	for (int i = 0; i < height * width * color; ++i) 
		data[i] += val;
}

template <class T>
void Image<T>::multiply(T val) {
	for (int i = 0; i < height * width * color; ++i) 
		data[i] *= val;
}

template <class T>
void Image<T>::set(T val) {
    for (int i = 0; i < height * width * color; ++i)
        data[i] = val;
}

template <class T>
void Image<T>::log() {
	for (int i = 0; i < height * width * color; ++i)
		data[i] = std::log(data[i]);
}

template <class T>
void Image<T>::exp() {
	for (int i = 0; i < height * width * color; ++i)
		data[i] = std::exp(data[i]);
}

template <class T>
void Image<T>::tonemap(Float key) {
    // Implement the simple equation 3 in Reinhard's SIGGRAPH '02 paper. 
    if (typeid(T) != typeid(Float)) {
        Log::warn() << "Failed to tone map. Only Float data type is supported." << endn;
        return;
    }
    if (color != 3) {
        Log::warn() << "Failed to tone map. Only three-color image is supported." << endn;
        return;
    }

    Float *luminance = new Float[height * width];

    // discard all NaNs if any    
    for (int i = 0; i < height * width * 3; ++i) {
        if (data[i] != data[i]) {
            data[i] = 0;            
        }
    }

    for (int i = 0; i < height * width; ++i) {
        luminance[i] = 0.27 * data[3 * i] + 0.67 * data[3 * i + 1] + 0.06 * data[3 * i + 2];                
    }
    
    // find log-average luminance
    Float delta = 1e-4; // NOTE: change this parameter to control the saturation of the tone mapped image.
    Float avg_luminance = 0.;    
    for (int i = 0; i < height * width; ++i) {        
        Float log = std::log(delta + luminance[i]);
        avg_luminance += log;        
    }

    avg_luminance = std::exp(avg_luminance / (height * width)); // equation 1

    Float scale = key / avg_luminance;

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
        data[3 * i    ] *= luminance[i];
        data[3 * i + 1] *= luminance[i];
        data[3 * i + 2] *= luminance[i];
    }
    delete [] luminance;
}

template <class T>
void Image<T>::tonemap(int x0, int y0, int x1, int y1, Float key) {
    // TODO: handle x0, x1, y0, y1
    // Implement the simple equation 3 in Reinhard's SIGGRAPH '02 paper. 
    if (typeid(T) != typeid(Float)) {
        Log::warn() << "Failed to tone map. Only Float data type is supported." << endn;
        return;
    }
    if (color != 3) {
        Log::warn() << "Failed to tone map. Only three-color image is supported." << endn;
        return;
    }

    Float *luminance = new Float[height * width];

    // discard all NaNs if any    
    for (int i = 0; i < height * width * 3; ++i) {
        if (data[i] != data[i]) {
            data[i] = 0;            
        }
    }

    for (int i = 0; i < height * width; ++i) {
        luminance[i] = 0.27 * data[3 * i] + 0.67 * data[3 * i + 1] + 0.06 * data[3 * i + 2];                
    }
    
    // find log-average luminance
    Float delta = 1e-4; // NOTE: change this parameter to control the saturation of the tone mapped image.
    Float avg_luminance = 0.;    
    for (int i = 0; i < height * width; ++i) {        
        Float log = std::log(delta + luminance[i]);
        avg_luminance += log;        
    }

    avg_luminance = std::exp(avg_luminance / (height * width)); // equation 1

    Float scale = key / avg_luminance;

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
        data[3 * i    ] *= luminance[i];
        data[3 * i + 1] *= luminance[i];
        data[3 * i + 2] *= luminance[i];
    }
    delete [] luminance;
}


static vector<Vec2> bressenham(const Vec2 &from, const Vec2 &to) {
    vector<Vec2> points;
    #define SWAP(x, y) (x ^= y ^= x ^= y)

    int x0 = from.x();
    int y0 = from.y();
    int x1 = to.x();
    int y1 = to.y();

    int Dx = x1 - x0; 
    int Dy = y1 - y0;
    int steep = (std::abs(Dy) >= std::abs(Dx));
    if (steep) {
        SWAP(x0, y0);
        SWAP(x1, y1);
        // recompute Dx, Dy after swap
        Dx = x1 - x0;
        Dy = y1 - y0;
    }
    int xstep = 1;
    if (Dx < 0) {
        xstep = -1;
        Dx = -Dx;
    }
    int ystep = 1;
    if (Dy < 0) {
        ystep = -1;              
        Dy = -Dy; 
    }
    int TwoDy = 2*Dy; 
    int TwoDyTwoDx = TwoDy - 2*Dx; // 2*Dy - 2*Dx
    int E = TwoDy - Dx; //2*Dy - Dx
    int y = y0;
    int xDraw, yDraw;    
    for (int x = x0; x != x1; x += xstep) {         
        if (steep) {                     
            xDraw = y;
            yDraw = x;
        } else {                 
            xDraw = x;
            yDraw = y;
        }       
        points.push_back(Vec2(xDraw, yDraw));

        // next
        if (E > 0) {
            E += TwoDyTwoDx; //E += 2*Dy - 2*Dx;
            y = y + ystep;
        } else {
            E += TwoDy; //E += 2*Dy;
        }
    }    
    return points;
}

template <class T>
void Image<T>::draw_line(const Vec2 &from, const Vec2 &to, const Rgb &color, Float alpha) {
    vector<Vec2> points = bressenham(from, to);

    for (int k = 0; k < points.size(); ++k) {
        int i = points[k].y() * width + points[k].x();        

        data[3 * i    ] = color.red()   * alpha + data[3 * i    ] * (1 - alpha);
        data[3 * i + 1] = color.green() * alpha + data[3 * i + 1] * (1 - alpha);
        data[3 * i + 2] = color.blue()  * alpha + data[3 * i + 2] * (1 - alpha);
    }
}

template <class T>
void Image<T>::set(const Vec2 &pixel, const Rgb &val) {
	int i = pixel.y();
	int j = pixel.x();	

	(*this)(i, j, 0) = val.red();
	(*this)(i, j, 1) = val.green();
	(*this)(i, j, 2) = val.blue();
}

template <class T>
void Image<T>::accumulate(const Vec2 &pixel, const Rgb &val) {
	int i = pixel.y();
	int j = pixel.x();	

	(*this)(i, j, 0) += val.red();
	(*this)(i, j, 1) += val.green();
	(*this)(i, j, 2) += val.blue();
}

template <class T>
bool Image<T>::is_valid_pixel(const Vec2 &pixel) const {
	int i = pixel.y();
	int j = pixel.x();
	if (i < 0 || i >= height || j < 0 || j >= width) return false;
	return true;
}

template <class T>
void Image<T>::flip_vertical() {
    T *new_data = new T[height * width * color];

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < color; ++k) {
                int source = color * (i                * width + j) + k;
                int dest   = color * ((height - i - 1) * width + j) + k;

                new_data[dest] = data[source];
            }
        }
    }

    delete [] data;
    data = new_data;
}

template <class T>
Rgb Image<T>::lookup(Float u, Float v) const {
    // nearest neighbor
    int i = (int)(v * (height - 1));
    int j = (int)(u * (width - 1));
    if (color == 1) {
        return Rgb(data[i * width + j]);
    } else if (color == 3) {
        return Rgb(data[3 * (i * width + j)    ], 
                   data[3 * (i * width + j) + 1],
                   data[3 * (i * width + j) + 2]);
    } else {
        return DefaultRgb::black;
    }
}

} // end namespace Renzoku
#endif
