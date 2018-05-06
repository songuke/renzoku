#ifndef _FILM_H_
#define _FILM_H_

#include "image.h"

namespace Renzoku {

class Film {

public:

    Film() {
    }

    Film(Size2 size) {
        img.initialize(0.0f, size.height, size.width, 3);
    }

    void initialize(Size2 size) {
        img.initialize(0.0f, size.height, size.width, 3);
    }

    ImageFloat &get_image() { return img; }

    void splat(Vec2 pixel, Rgb color) {
        int i = (int)pixel.y();
        int j = (int)pixel.x();
        img(i, j, 0) += color.red();
        img(i, j, 1) += color.green();
        img(i, j, 2) += color.blue();
    }

protected:
    ImageFloat img;

    // TODO: implement reconstruction filter
};

}   // end namespace Renzoku

#endif