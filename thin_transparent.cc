#include "thin_transparent.h"
#include "onb.h"

namespace Renzoku {

ThinTransparent::ThinTransparent(const Rgb &color) : color(color) {
}

Rgb ThinTransparent::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    wi = -wo;
    pdf = 1.0f;
    return color * INV_PI / fabs(dot(wi, uvn.n()));
}

Rgb ThinTransparent::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf, Rgb &reflectivity) {
    wi = -wo;
    pdf = 1.0f;
    return color * INV_PI / fabs(dot(wi, uvn.n()));;
}

Rgb ThinTransparent::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    return DefaultRgb::black;
}

Float ThinTransparent::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    return 0.0f;
}


} // end namespace Renzoku