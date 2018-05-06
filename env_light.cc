#include "env_light.h"
#include "surface.h"
#include "sphere.h"
#include "image.h"

namespace Renzoku {

EnvLight::EnvLight(Surface *surface, ImageFloat *envmap) 
    : surface(surface), envmap(envmap) {
    sphere = (Sphere *)surface->get_shape();

    // align to the world's z axis. 
    uvn = new Onb();
    uvn->set_identity();

    // build a discrete pdf
    pdf2d.set_distribution(*envmap);
}

EnvLight::~EnvLight() {
    if (uvn)
        delete uvn;
}

Float EnvLight::area() const {
    return sphere->area();
}

Rgb EnvLight::power() const {
    return Rgb(sphere->area(), sphere->area(), sphere->area());
}
    
bool EnvLight::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    return false;
}

bool EnvLight::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return false;
}

void EnvLight::sample(Random &rd, Vec3 &wo, Float &pdf) const {
    /*
    // uniformly sample an out-going direction from the center of the hemisphere
    Vec2 uv;
    uv.random(rd);
    Float phi = 2 * A_PI * uv[0];
    Float cos_theta = 1.0 - 2 * uv[1];
    Float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    
    // up vector is Z
    wo = Vec3(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta);

    // rotate to world
    wo = uvn->local_to_world(wo);
    
    pdf = INV_4PI;
    */

    // importance sampling based on intensity in the environment map
    Float u, v, pdf_uv;
    int row, col;
    pdf2d.sample(rd, row, col, pdf_uv);

    Log::info() << row << " " << col << endn;

    v = (Float)row / envmap->get_height();
    u = (Float)col / envmap->get_width();

    Float phi = TWO_PI * u;
    Float theta = A_PI * v;

    Float sin_theta = sin(theta);
    wo = Vec3(sin_theta * cos(phi), sin_theta * sin(phi), cos(theta));

    // rotate to world
    wo = uvn->local_to_world(wo);
    pdf = pdf_uv / (INV_2PI * INV_PI * sin_theta);
}

Float EnvLight::pdf(const Vec3 &wo) const {
    Vec3 d = uvn->world_to_local(wo);
    
    Float theta = acos(d.z());
    Float phi = atan2(d.y(), d.x());
    
    Float u = (phi + A_PI) * INV_2PI;
    Float v = theta * INV_PI;

    int row = v * envmap->get_height();
    int col = u * envmap->get_width();

    return pdf2d.probability(row, col);
}

void EnvLight::sample(Random &rd, Vec3 &wo, Rgb &radiance, Float &pdf) const {
    this->sample(rd, wo, pdf);

    this->query(wo, radiance);
}


void EnvLight::sample(Random &rd, const Vec3 &n, Vec3 &wo, Float &pdf) const {
    // cosine weighted sampling
    Vec2 sampler;
    sampler.random(rd);
        
    Float phi = 2 * A_PI * sampler.x();
    Float sin_theta = sqrt(sampler.y());
    Float x = sin_theta * cos(phi);
    Float y = sin_theta * sin(phi);
    Float z = sqrt(1 - sin_theta * sin_theta);
 
    Onb basis;
    basis.init_from_n(n);

    wo = unit_vector(basis.local_to_world(Vec3(x, y, z)));
    
    pdf = z * INV_PI;
}
void EnvLight::pdf(const Vec3 &n, const Vec3 &wo, Float &pdf) const {
    pdf = fabs(dot(n, wo)) * INV_PI;
}

void EnvLight::sample(Random &rd, Vec3 &wo, Float &pdf, const Receiver &patch) const {
    this->sample(rd, wo, pdf, patch);
}

void EnvLight::sample(Random &rd, Vec3 &wo, Rgb &radiance, Float &pdf, const Receiver &patch) const {
    this->sample(rd, wo, radiance, pdf, patch);
}

Float EnvLight::pdf(const Vec3 &wo, const Receiver &patch) const {
    return this->pdf(wo, patch);
}

void EnvLight::query(const Vec3 &wo, Rgb &radiance) const {
    // texture rotation
    Vec3 d = uvn->world_to_local(wo);

    Float theta = acos(d.z());
    Float phi = atan2(d.y(), d.x());
    
    Float u = (phi + A_PI) * INV_2PI;
    Float v = theta * INV_PI;
    
    radiance = envmap->lookup(u, v);
}

Rgb EnvLight::first_bounce(Scene *scene, const Receiver &r, DirectionSampler *dir_sampler) const {
    return DefaultRgb::black;
}

Rgb EnvLight::irradiance(Scene *scene, const LocalGeometry &dg, LightSamplingRecord &sr) const {
    sr.is_valid = false;
    return DefaultRgb::black;
}

} // end namespace Renzoku
