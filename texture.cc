#include "texture.h"

#include "image.h"

namespace Renzoku {

Texture::Texture(const string img) {
    buf = new ImageFloat(img.c_str());
}

Texture::Texture(ImageFloat &img) {
    buf = &img;
}

Texture::~Texture() {
    if (buf)
        delete buf;
}

Rgb Texture::lookup(const Vec2 &uv) {
    // simple repeat mode
    Float col = uv.x() - (int)uv.x();
    Float row = uv.y() - (int)uv.y();

    col = fabs(col);
    row = fabs(row);
    return buf->lookup(col, row);
}

} // end namespace