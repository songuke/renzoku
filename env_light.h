#ifndef _ENV_LIGHT_H_
#define _ENV_LIGHT_H_

#include "common.h"
#include "light.h"
#include "pdf.h"

namespace Renzoku {

class EnvLight : public Light {
public:
    EnvLight(Surface *_surface, ImageFloat *envmap);
    virtual ~EnvLight();

    virtual Light::Type get_light_type() const {
        return Light::ENV_LIGHT;
    }

    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

    virtual Float area() const;
    virtual Rgb power() const;
    
    /** 
     * Uniform sample the sphere
     */
    virtual void sample(Random &rd, Vec3 &wo, Float &pdf) const;        
    virtual void sample(Random &rd, Vec3 &wo, Rgb &radiance, Float &pdf) const;
    
    /**
     * Perform the cosine weighted sampling of the hemisphere the orientation n belongs to.
     */
    virtual void sample(Random &rd, const Vec3 &n, Vec3 &wo, Float &pdf) const;
    virtual void pdf(const Vec3 &n, const Vec3 &wo, Float &pdf) const;

    virtual void sample(Random &rd, Vec3 &wo, Float &pdf, const Receiver &patch) const;
    virtual void sample(Random &rd, Vec3 &wo, Rgb &radiance, Float &pdf, const Receiver &patch) const;
    
    virtual Float pdf(const Vec3 &wo) const;            
    virtual Float pdf(const Vec3 &wo, const Receiver &patch) const;

    virtual void query(const Vec3 &wo, Rgb &radiance) const;

    virtual Rgb first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler = NULL) const;
    virtual Rgb irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const;

protected:
    ImageFloat *envmap;
    Surface *surface;
    Sphere *sphere;
    Onb *uvn;       // local transformation (rotation) of the sphere

    DiscretePdf2D pdf2d;
};

} // end namespace Renzoku

#endif
