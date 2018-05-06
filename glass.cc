#include "glass.h"
#include "onb.h"
#include "random.h"
#include "log.h"

namespace Renzoku {

Glass::Glass(Float inner) : inner(inner), outer(1.f), kd(DefaultRgb::white) {
}
    
Glass::Glass(Float inner, Float outer) : inner(inner), outer(outer), kd(DefaultRgb::white) {
}
   
Glass::Glass(Float inner, const Rgb &kd) : inner(inner), outer(1.f), kd(kd) {
}

Glass::Glass(Float inner, Float outer, const Rgb &kd) : inner(inner), outer(outer), kd(kd) {
}

/**
 * Convention:
 *
 * Glass normal vector is always facing to outer environment.
 *
 * Glass surface is not two-sided. In other words, glass objects should be modelled as thick object. 
 */
Rgb Glass::sample(Random &rd, const Vec3 &nn, Float nr, Float nt, const Vec3 &wo, Vec3 &wi, Float &pdf) {

    Float dot_wn = dot(wo, nn);
    
    Float n1 = nr;
    Float n2 = nt;

    if (dot_wn < 0.0f) {
        n1 = nt;
        n2 = nr;
    }
                
    Float n1n2 = n1 / n2;       
    Float delta = 1.f - (n1n2 * n1n2) * (1.f - dot_wn * dot_wn);
    if (delta <= ZERO_EPSILON) {                                      // total internal reflection
        
        wi = 2.0f * dot_wn * nn - wo;
        pdf = 1.f;
                
        if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
            Float f = 1.0f / fabs(dot(wi, nn));                       // reflection BSDF based on delta distribution
            return kd * f * INV_PI;
        } else return DefaultRgb::black;

    } else { 
        
        // compute Fresnel term to determine the amount of energy that reflects
        // 
        // total energy:
        // E = F * reflection + (1 - F) * transmission
        // 
        // Each time we choose either reflection or transmission to trace
        // with probability F for reflection and 1 - F for transmission.
        // 
        // The estimator of E is:
        // E' = F * reflection / pdf            with probability F (so pdf = F and F cancels)
        // E' = (1 - F) * transmission / pdf    with probability 1 - F
        //
        // The average value of estimator E' is F * reflection + (1 - F) * transmission, which is the value we want.
        //

        wi = 2.0f * dot_wn * nn - wo;

        Float F0 = (n1 - n2) / (n1 + n2);
        F0 *= F0;
            
        // due to perfect reflection, the half vector coincides with the normal
        // so, F = F0.
        /*
        Vec3 wh = unit_vector(wi + wo);
        Float p = 1.f - dot(nn, wh);
        Float F = F0 + (1.f - F0) * p * p * p * p * p;
        */
        Float F = F0;
                
        Float u = rd();                                             // follow either reflection or transmission 
        if (u < F) {
            
            pdf = F;
            if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
                Float f = 1.0f / fabs(dot(wi, nn));                
                return F * kd * f * INV_PI;
            } else {                
                return DefaultRgb::black;
            }

        } else { 
            
            if (dot_wn > 0.0f) {
                wi = - n1n2 * (wo - dot_wn * nn) - sqrt(delta) * nn;       
            } else {
                wi = - n1n2 * (wo - dot_wn * nn) + sqrt(delta) * nn;
            }

            if (dot(wi, wo) > 0.0f) {
                pdf = 0.0f;
                wi = Vec3(0.0f, 0.0f, 0.0f);
                
                Log::info() << "Glass: wi, wo same direction" << endn;
                return DefaultRgb::black;                
            }
                        
            pdf = 1.f - F;

            if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
                Float f = (n2*n2 / n1*n1) / fabs(dot(wi, nn));        // refraction BSDF based on delta distribution
                return (1.f - F) * kd * f * INV_PI;
            } else {
                return DefaultRgb::black;
            }

        }
    }    
}

Rgb Glass::sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) {
    return sample(rd, uvn.n(), outer, inner, wo, wi, pdf);
}

Rgb Glass::eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    
    Vec3 nn = uvn.n();
    Float dot_wn = dot(wo, nn);
    
    Float n1 = outer;
    Float n2 = inner;

    if (dot_wn < 0.0f) {
        n1 = inner;
        n2 = outer;
    }
    
    Float n1n2 = n1 / n2;       

    Float delta = 1.f - (n1n2 * n1n2) * (1.f - dot_wn * dot_wn);
    if (delta <= ZERO_EPSILON) {

        Vec3 mirror_wi = 2.0f * dot_wn * nn - wo;
        if (1.0f - fabs(dot(wi, mirror_wi)) <= ZERO_EPSILON) {
            if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
                Float f = 1.0f / fabs(dot(wi, nn));                 
                return kd * f * INV_PI;   
            }
        }

    } else {

        Float F0 = (n1 - n2) / (n1 + n2);
        F0 *= F0;            
        Float F = F0;

        Vec3 mirror_wi = 2.0f * dot_wn * nn - wo;
        if (1.0f - fabs(dot(wi, mirror_wi)) <= ZERO_EPSILON) {
            if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
                Float f = 1.0f / fabs(dot(wi, nn));                
                return F * kd * f * INV_PI;
            }
        }
        
        Vec3 transmit_wi;
        if (dot_wn >= 0.0f) 
            transmit_wi = - n1n2 * (wo - dot_wn * nn) - sqrt(delta) * nn;
        else
            transmit_wi = - n1n2 * (wo - dot_wn * nn) + sqrt(delta) * nn;

        if (1.0f - fabs(dot(wi, transmit_wi)) <= ZERO_EPSILON) {
            if (fabs(dot(wi, nn)) > ZERO_EPSILON) {
                Float f = (n2*n2 / n1*n1) / fabs(dot(wi, nn));                
                return (1.f - F) * kd * f * INV_PI;
            }
        }

    }

    return DefaultRgb::black;
}

/** 
 * This pdf cannot just return 0.0f since we need to pdf(wi) "as if"
 * wi is sampled from wo. 
 *
 * This is important for MIS weight calculation in bidir path tracing.
 */
Float Glass::pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) {
    
    Vec3 nn = uvn.n();
    Float dot_wn = dot(wo, nn);
    
    Float n1 = outer;
    Float n2 = inner;

    if (dot_wn < 0.0f) {
        n1 = inner;
        n2 = outer;
    }
    
    Float n1n2 = n1 / n2;       

    Float delta = 1.f - (n1n2 * n1n2) * (1.f - dot_wn * dot_wn);
    if (delta <= ZERO_EPSILON) {

        Vec3 mirror_wi = 2.0f * dot_wn * nn - wo;
        if (1.0f - fabs(dot(wi, mirror_wi)) <= ZERO_EPSILON) return 1.0f;

    } else {

        Float F0 = (n1 - n2) / (n1 + n2);
        F0 *= F0;            
        Float F = F0;
        
        Vec3 mirror_wi = 2.0f * dot_wn * nn - wo;
        if (1.0f - fabs(dot(wi, mirror_wi)) <= ZERO_EPSILON) return F;
        
        Vec3 transmit_wi;
        if (dot_wn >= 0.0f) 
            transmit_wi = - n1n2 * (wo - dot_wn * nn) - sqrt(delta) * nn;
        else
            transmit_wi = - n1n2 * (wo - dot_wn * nn) + sqrt(delta) * nn;

        if (1.0f - fabs(dot(wi, transmit_wi)) <= ZERO_EPSILON) return 1.0f - F;

    }
    return 0.0f;
}

} // end namspace Renzoku
