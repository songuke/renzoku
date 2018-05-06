#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "common.h"

namespace Renzoku {

class Texture {
public:
    Texture(const string img);
    Texture(ImageFloat &img);
    virtual ~Texture();

    Rgb lookup(const Vec2 &uv);

    inline ImageFloat *get_buffer() const;

protected:
    ImageFloat *buf;    
};

inline ImageFloat* Texture::get_buffer() const {
    return buf;
}

} // end namespace Renzoku

#endif
