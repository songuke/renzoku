#include "mirror.h"
#include "onb.h"
#include "random.h"

namespace Renzoku {
Mirror::Mirror() : kd(DefaultRgb::white) {         
}

Mirror::Mirror(const Rgb &kd) : kd(kd) {    
}

/**
 * Convention:
 *
 * in, out are out-going vectors from the intersection point.
 */
Rgb Mirror::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    Vec3 nn = uvn.n();
    Float dot_wn = dot(wo, nn);
    
    if (dot_wn <= ZERO_EPSILON) {
        wi = Vec3(0.0f, 0.0f, 0.0f);
        pdf = 0.0f;
        return DefaultRgb::black;
    }    
    
    wi = 2.0f * dot_wn * nn - wo;     // wi is still a unit vector
    pdf = 1.0f;
    return kd / dot(wi, nn);
}

Rgb Mirror::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    Vec3 nn = uvn.n();
    Float dot_wn = dot(wo, nn);
    if (dot_wn <= ZERO_EPSILON)
        return DefaultRgb::black;

    Vec3 wr = 2.0f * dot_wn * nn - wo;
    if (fabs(dot(wr, wi)) >= 1.0f - ZERO_EPSILON) {
        // reflection BSDF based on delta distribution
        return kd / dot(wi, uvn.n());
    }
    return DefaultRgb::black;    
}

Float Mirror::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    Vec3 nn = uvn.n();
    Float dot_wn = dot(wo, nn);
    if (dot_wn <= ZERO_EPSILON)
        return 0.0f;

    Vec3 wr = 2.0f * dot_wn * nn - wo;
    if (fabs(dot(wr, wi)) >= 1.0f - ZERO_EPSILON) {
        // reflection BSDF based on delta distribution
        return 1.0f;
    }
    return 0.0f;
}

}
