#include "photon_light.h"

namespace Renzoku {
PhotonLight::PhotonLight(Photon *photon, Material *m)
    : photon(photon), material(m), radius(0.)
{    
}

void PhotonLight::set_radius(Float radius) {
    this->radius = radius;
}

bool PhotonLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return false;
}

bool PhotonLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return false;
}

Float PhotonLight::area() const {
    return 0.;
}

Rgb PhotonLight::power() const {
    return photon->flux;
}

void PhotonLight::sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const {
    throw "Not implemented for photon light.";
}

void PhotonLight::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    throw "Not implemented for photon light.";
}

void PhotonLight::query(const Vec3 &p, const Vec3 &wo, Rgb &radiance) {
    throw "Not implemented for photon light.";
}

Float PhotonLight::pdf(const Vec3 &p) const {
    return 0.0f;
}

Float PhotonLight::pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo) const {
    return 0.0f;
}

 
} // end namespace Renzoku
