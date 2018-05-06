#include "phong.h"
#include "log.h"

namespace Renzoku {

ModifiedPhong::ModifiedPhong(Rgb _kd, Rgb _ks, Float _n) : kd(_kd), ks(_ks), n(_n) {     
    rho_d = kd.max();
    rho_s = ks.max();
    
    // normalize rho_d and rho_s for sampling
    Float total = rho_d + rho_s;
	if (total > 0.0f) {
		rho_d /= total;
		rho_s /= total;
	}
}

/**
 * Only sample a ray in a component (diffuse, specular, or none)
 */
Rgb ModifiedPhong::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    if (dot(wo, uvn.n()) <= ZERO_EPSILON) {
        wi = Vec3(0.0f, 0.0f, 0.0f);
        pdf = 0.0f;
        return DefaultRgb::black;
    }

    Float epsilon = rd();
        
    if (epsilon < rho_d) {
        Vec2 sampler;
        sampler.random(rd);

        Float phi = TWO_PI * sampler.x();
        Float sin_theta = sqrt(sampler.y());
    
        Float x = sin_theta * cos(phi);
        Float y = sin_theta * sin(phi);
        Float z = sqrt(1.0f - sin_theta * sin_theta);
 
        wi = unit_vector(uvn.local_to_world(Vec3(x, y, z)));

    } else {

        // sample pdf = (n + 1) / 2PI * cos(alpha)^n
        Vec2 sampler;
        sampler.random(rd);

        // choose an easy coordinate system
        Onb onb;
        onb.init_from_nu(wo, uvn.n());
        
        Float phi = TWO_PI * sampler.x();
        Float alpha = acos(pow(sampler.y(), 1.0f / (n + 1)));

        Vec3 wr = onb.spherical_to_world(alpha, phi);

        Vec3 nn = uvn.n();
        Float dot_wr_nn = dot(wr, nn);        
        if (dot_wr_nn < 0.0f) {
            // reject this sample since it will generate wi under surface
            wi = Vec3(0.0f, 0.0f, 0.0f);
            pdf = 0.0f;
            return DefaultRgb::black;
        }
        wi = 2 * dot_wr_nn * nn - wr;
    }

    pdf = this->pdf(uvn, wo, wi);
    Rgb brdf = this->eval(uvn, wo, wi);
    return brdf;
}

Rgb ModifiedPhong::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    if (dot(wi, uvn.n()) <= ZERO_EPSILON || dot(wo, uvn.n()) <= ZERO_EPSILON)
        return DefaultRgb::black;

    Rgb brdf = kd * INV_PI;

    Vec3 nn = uvn.n();
    Vec3 wr = 2 * dot(wi, nn) * nn - wi;
    Float cos_alpha = dot(wo, wr);
    if (cos_alpha > 0.0f) {
        brdf += ks * (n + 2) * INV_2PI * pow(cos_alpha, n);        
    }
    return brdf;
}

Float ModifiedPhong::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    if (dot(wi, uvn.n()) <= ZERO_EPSILON || dot(wo, uvn.n()) <= ZERO_EPSILON)
        return 0.0f;
    
    // the "average" pdf becayse a ray could be sampled by both parts
    Vec3 nn = uvn.n();
    Float pdf = rho_d * (dot(nn, wi) * INV_PI);

    Vec3 wr = 2 * dot(wi, nn) * nn - wi;
    Float cos_alpha = dot(wo, wr);
    if (cos_alpha > 0.0f) {
        pdf += rho_s * (n + 1) * INV_2PI * pow(cos_alpha, n);
    }
    return pdf;    
}

}  // end namespace
