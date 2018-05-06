#include "ward.h"
#include "log.h"

namespace Renzoku {

Ward::Ward(Rgb _kd, Rgb _ks, Float _ax, Float _ay) : kd(_kd), ks(_ks), ax(_ax), ay(_ay) { 
    ax2 = ax * ax;
    ay2 = ay * ay;
    ay_ax = ay / ax;
    inv_four_pi_ax_ay = 1.0f / (4 * A_PI * ax * ay);

    // for simplicity, let rho_s = ks (maximum possible value - which may result in oversampling
    // in specular term).
    // to be more accurate, need to sample an outgoing direction and estimate rho_s. 
    rho_d = kd.max();
    rho_s = ks.max();
    
    // normalize rho_d and rho_s for sampling
    Float total = rho_d + rho_s;
	if (total > 0.0f) {
		rho_d /= total;
		rho_s /= total;
	}
}

Rgb Ward::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
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
        Float z = sqrt(1 - sin_theta * sin_theta);
 
        wi = unit_vector(uvn.local_to_world(Vec3(x, y, z)));
        
    } else {
        // requires u, v in (0, 1) (zero should not be included)
        Vec2 sampler(rd.rand_exclude_01(), rd.rand_exclude_01());
        
        Vec3 local_wo = unit_vector(uvn.world_to_local(wo));
                
        Float angle = TWO_PI * sampler.x();
        Float phi = atan(ay_ax * tan(angle));

        // atan() returns [-HALF_PI, HALF_PI] while phi needs to be in [0, TWO_PI].        
        // We need to ensure phi ranges in [0, TWO_PI] and to be in the same quadrant with TWO_PI * sampler.x()
        // If this is not considered, specular highlight that falls on triangle boundaries splits into two halves. 
        // Disable the code below and use Veach's MIS scene to see. 
        /*
        if (-HALF_PI <= phi && phi <= 0) {
            if (HALF_PI <= angle && angle <= A_PI) { // opposite quadrant
                phi += A_PI;
            } else if (A_PI + HALF_PI <= angle && angle <= TWO_PI) {
                // same quadrant
            } else {
                Log::error() << "Ward BRDF: angle in error quadrant 1" << endl;
            }
        } else if (0 <= phi && phi <= HALF_PI) {
            if (A_PI <= angle && angle <= A_PI + HALF_PI) { // opposite quadrant
                phi += A_PI;
            } else if (0 <= angle && angle <= HALF_PI) {
                // same quadrant
            } else {
                Log::error() << "Ward BRDF: angle in error quadrant 2" << endl;
            }
        } else {
            Log::error() << "Ward BRDF: phi is not -pi/2, pi/2: " << phi << endl;
        }*/

        if (sin(phi) * sin(angle) < 0 || cos(phi) * cos(angle) < 0)
            phi += A_PI;

        Float cos_phi = cos(phi);
        Float sin_phi = sin(phi);
        Float theta = atan(sqrt(-log(sampler.y()) / 
                      (cos_phi * cos_phi / ax2 + sin_phi * sin_phi / ay2) ));
        
        Float sin_theta = sin(theta);
        Float cos_theta = cos(theta);
        Vec3 wh(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
                
        Vec3 local_wi = 2 * dot(local_wo, wh) * wh - local_wo;        
        wi = uvn.local_to_world(local_wi);
    }
                                            
    pdf = this->pdf(uvn, wo, wi);           // probability of wi must account for both components
    Rgb brdf = this->eval(uvn, wo, wi);     // brdf is the kd + ks
                                            // the rho_d, rho_s are for sampling wi, not for computing brdf
    return brdf;
}

Rgb Ward::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    Vec3 local_wi = uvn.world_to_local(wi);
    Vec3 local_wo = uvn.world_to_local(wo);

    if (local_wi.z() <= ZERO_EPSILON || local_wo.z() <= ZERO_EPSILON)
        return DefaultRgb::black;
        
    Vec3 wh = local_wi + local_wo;
    Float expo = exp(-(wh.x() * wh.x() / ax2 + wh.y() * wh.y() / ay2) / (wh.z() * wh.z()));
    Rgb brdf = kd * INV_PI + ks * inv_four_pi_ax_ay / sqrt(local_wi.z() * local_wo.z()) * expo;

    return brdf;    
}

Float Ward::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    if (dot(wi, uvn.n()) <= ZERO_EPSILON || dot(wo, uvn.n()) <= ZERO_EPSILON)
        return 0.0f;
        
    Float z = dot(uvn.n(), wi);
        
    Vec3 wh = unit_vector(uvn.world_to_local(wi + wo));
    Float expo = exp(-(wh.x() * wh.x() / ax2 + wh.y() * wh.y() / ay2) / (wh.z() * wh.z()));
    Float cos_theta = wh.z();
    
    Vec3 local_wo = unit_vector(uvn.world_to_local(wo));
    Float cos_theta3 = cos_theta * cos_theta * cos_theta;
        
    // average probability
    return rho_d * (z * INV_PI) +
           rho_s * (inv_four_pi_ax_ay / (dot(wh, local_wo) * cos_theta3) * expo);
}

} // end namespace Renzoku
