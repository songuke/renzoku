#include "lambertian.h"

namespace Renzoku {

Rgb Lambertian::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    if (dot(wo, uvn.n()) <= ZERO_EPSILON) {
        wi = Vec3(0.0f, 0.0f, 0.0f);
        pdf = 0.0f;
        return DefaultRgb::black;
    }

    Vec2 sampler;
    sampler.random(rd);

    // diffuse BRDF is independent of incident ray
    Float phi = TWO_PI * sampler.x();
    Float sin_theta = sqrt(sampler.y());
    Float x = sin_theta * cos(phi);
    Float y = sin_theta * sin(phi);
    Float z = sqrt(1 - sin_theta * sin_theta);
 
    wi = unit_vector(uvn.local_to_world(Vec3(x, y, z)));
    
    pdf = z * INV_PI; // cos_theta / PI
    return kd * INV_PI;
}

/*
Rgb Lambertian::sample_uniform(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    if (dot(wo, uvn.n()) <= ZERO_EPSILON) {
        wi = Vec3(0.0f, 0.0f, 0.0f);
        pdf = 0.0f;
        return DefaultRgb::black;
    }

    Vec2 sampler;
    sampler.random(rd);
            
    Float phi = 2 * A_PI * sampler.x();
    Float cos_theta = sampler.y();
    Float sin_theta = sqrt(1. - cos_theta * cos_theta); // PBRT adds a max(0, 1 - cos_theta^2). Why?

    Float x = sin_theta * cos(phi);
    Float y = sin_theta * sin(phi);
    Float z = sampler.y(); // cos(theta)
 
    wi = unit_vector(uvn.local_to_world(Vec3(x, y, z)));

    // due to numerical accuracy, wi still can be under hemisphere, so need to double check before return
    if (dot(wi, uvn.n()) < 0) {
        cout << "Before: " << z << endl;
        cout << "After:  " << dot(wi, uvn.n()) << endl;
        //cout << "ERROR: diffuse sampling vector under hemisphere." << endl;
        wi = Vec3(0, 1, 0);
        pdf = 1;
        return Rgb(0, 0, 0);
    }
    pdf = 1.0f / TWO_PI;
    
    return kd * INV_PI;
}*/

Rgb Lambertian::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    if (dot(wi, uvn.n()) <= ZERO_EPSILON || dot(wo, uvn.n()) <= ZERO_EPSILON)
        return DefaultRgb::black;

    return kd * INV_PI;
}

Float Lambertian::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    if (dot(wi, uvn.n()) <= ZERO_EPSILON || dot(wo, uvn.n()) <= ZERO_EPSILON)
        return 0.0f;

    return dot(uvn.n(), wi) * INV_PI;
}

} // end namespace Renzoku
