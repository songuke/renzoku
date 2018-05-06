#ifndef _PHOTON_LIGHT_
#define _PHOTON_LIGHT_

#include "common.h"
#include "light.h"
#include "photon.h"

namespace Renzoku {
class PhotonLight : public Light {
public:
    PhotonLight(Photon *photon, Material *m);

    virtual Light::Type get_light_type() const {
        return Light::PHOTON_LIGHT;
    }

    bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    Float area() const;
    Rgb power() const;
    void set_radius(Float radius);
    
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Vec3 &wo, Rgb &radiance, Float &pdf_p, Float &pdf_wo) const;
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;

    virtual Float pdf(const Vec3 &p) const;
    virtual Float pdf(const Vec3 &p, const Vec3 &n, const Vec3 &wo) const;

    void query(const Vec3 &p, const Vec3 &wo, Rgb &radiance);
    
public:
    Float center;
    Float radius;
    Photon *photon;
    Material *material;
};

typedef vector<PhotonLight *> PhotonLights;

}; // end namespace Renzoku

#endif

