#ifndef _VPL_INTERFACE_H_
#define _VPL_INTERFACE_H_

#include "brdf_point_light.h"

namespace Renzoku {

class IVplGenerator;

/**
 * This pseudo abstract class defines the common interface that 
 * each VPL implementation should support. 
 */ 
class IVirtualPointLight {
public:
    IVirtualPointLight() : vplgen(NULL), num_passes(1), cur_pass(0), max_vpls(0), num_vpls_save_interval(0) {
    }
    
    inline void set_start_particles(int num_particles);
    inline void set_vpl_generator(IVplGenerator *gen);
    inline void set_max_vpls(int max_vpls);
    inline void set_num_passes(int pass);
    inline void set_save_every(int vpls);

protected:
    bool is_store_direct_vpls;

    int num_start_particles;  
    BrdfPointLights all_vpls;

    IVplGenerator *vplgen;

    int num_passes, cur_pass;
    int max_vpls;
    int num_vpls_save_interval;
};

class IVirtualSphericalLight {
public:
    IVirtualSphericalLight() : radius_scale(1.0f) {
    }

    inline void set_radius_scale(Float radius_scale);

protected:
    Float radius_scale;
};


struct RadianceTerms {
    Rgb surface_brdf;
    Rgb light_brdf;
    Float form_factor;    
    Rgb power;
        
    Rgb radiance_no_clamp;
    Rgb radiance;
    bool is_clamped;
    Float w_clamp;          // the weight that is equivalent to clamping

    RadianceTerms(const Rgb &power, const Rgb &surface_brdf, const Rgb &light_brdf, Float form_factor) {
        this->radiance = DefaultRgb::black;
        this->radiance_no_clamp = DefaultRgb::black;
        this->is_clamped = false;
        this->w_clamp = 1.0f;

        this->power = power;
        this->surface_brdf = surface_brdf;
        this->light_brdf = light_brdf;
        this->form_factor = form_factor;
    }

    void eval_clamp_form_factor(Float threshold) {
        if (form_factor > threshold) {
            is_clamped = true;        
            w_clamp = threshold / form_factor;
            radiance = light_brdf * surface_brdf * threshold * power;
        } else {
            is_clamped = false;
            radiance = light_brdf * surface_brdf * form_factor * power;
        }    
        radiance_no_clamp = power * form_factor * surface_brdf * light_brdf;
    }

    void eval_clamp_reflectivity(Float threshold) {
        Rgb reflectivity = form_factor * surface_brdf * light_brdf;

        Float luminance = reflectivity.value();        
        if (luminance > threshold) {
            is_clamped = true;

            w_clamp = threshold / luminance;
            radiance = w_clamp * reflectivity * power;            

        } else {
            is_clamped = false;
            radiance = reflectivity * power;
        }
        radiance_no_clamp = power * form_factor * surface_brdf * light_brdf;
    }

    void eval_clamp_radiance(Float threshold) {        
        this->radiance_no_clamp = power * form_factor * surface_brdf * light_brdf;
                
        Float brightness = radiance_no_clamp.value();        
        if (brightness > threshold) {
            is_clamped = true;
            
            w_clamp = threshold / brightness;
            radiance = w_clamp * radiance_no_clamp;
        }
        else {
            is_clamped = false;
            radiance = radiance_no_clamp;
        }
    }

    void eval_no_clamp() {
        radiance_no_clamp = power * form_factor * surface_brdf * light_brdf;
        radiance = radiance_no_clamp;
        is_clamped = false;
    }
};


class IVirtualPointLightEvaluator {
public:
    virtual Rgb radiance(Scene *scene, const Receiver &r, const BrdfPointLight &light) {
        bool visible;
        return radiance(scene, r, light, visible);
    }

    virtual Rgb radiance(Scene *scene, const Receiver &r, const BrdfPointLight &light, bool &visible) = 0;

    /**
     * Incident radiance from the VPL light. No visibility is evaluated.
     */
    virtual Rgb incident_radiance(Scene *scene, const Receiver &r, const BrdfPointLight &light) = 0;
};

/**
 * Method to calculate radiance from all VPLs.
 */
struct VplRadianceEstimate {
    enum Type {
        NONE,               // accmulate every VPL one by one
        BIAS_ONLY,          // estimate bias caused by VPLs
        LIGHTCUTS,          
        MRCS, 
        LIGHTSLICE,
        POWER_SAMPLING,
        LIGHTSLICE_GUIDED_PATH_TRACE,
        VORBA_GUIDED_PATH_TRACE
    };
};

struct VplRadianceClamp {
    enum Type {
        NONE, 
        FORM_FACTOR,
        REFLECTIVITY,
        RADIANCE
    };
};

void IVirtualPointLight::set_start_particles(int num_particles) {
    num_start_particles = num_particles;
}

inline void IVirtualPointLight::set_vpl_generator(IVplGenerator *gen) {
    vplgen = gen;
}

inline void IVirtualPointLight::set_num_passes(int pass) {
    this->num_passes = pass;
}

inline void IVirtualPointLight::set_max_vpls(int max_vpls) {
    this->max_vpls = max_vpls;
}

inline void IVirtualSphericalLight::set_radius_scale(Float radius_scale) {
    this->radius_scale = radius_scale;
}

inline void IVirtualPointLight::set_save_every(int vpls) {
    this->num_vpls_save_interval = vpls;
}

} // end namespace 

#endif