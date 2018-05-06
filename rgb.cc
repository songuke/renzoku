#include "rgb.h"

namespace Renzoku {

const Rgb DefaultRgb::white(1.0f, 1.0f, 1.0f);
const Rgb DefaultRgb::black(0.0f, 0.0f, 0.0f);
const Rgb DefaultRgb::grey(0.5f, 0.5f, 0.5f);

const Rgb DefaultRgb::red(1.0f, 0.0f, 0.0f);
const Rgb DefaultRgb::green(0.0f, 1.0f, 0.0f);
const Rgb DefaultRgb::blue(0.0f, 0.0f, 1.0f);

const Rgb DefaultRgb::yellow(1.0f, 1.0f, 0.0f);
const Rgb DefaultRgb::cyan(0.0f, 1.0f, 1.0f);
const Rgb DefaultRgb::magenta(1.0f, 0.0f, 1.0f);

Rgb Rgb::from_hsv(Float hue, Float saturation, Float value) {
    Float r = 0.f, g = 0.f, b = 0.f;

	while (hue > 360) hue -= 360;
	if (saturation == 0) {//gray
	  r = g = b = value;
	} else {
        int h;
        float f, p, q, t, section;
        section = hue / 60.f;
		h = (int)section;
		f = section - h;
		p = value * (1 - saturation);
		q = value * (1 - f*saturation);
		t = value * (1 - (1 - f)*saturation);
		
		switch (h) {
		case 0:
			r = value;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = value;
			b = p;
			break;
		case 2:
			r = p;
			g = value;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = value;
			break;
		case 4:
			r = t;
			g = p;
			b = value;
			break;
		case 5:
			r = value;
			g = p;
			b = q;
			break;
		}
		
	}
    return Rgb(r, g, b);
}

} // end namespace