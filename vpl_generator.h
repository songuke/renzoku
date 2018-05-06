#ifndef _VPL_GENERATOR_H_
#define _VPL_GENERATOR_H_

#include "brdf_point_light.h"

namespace Renzoku {

class IVplGenerator {    
public:    
    IVplGenerator() : is_store_direct_vpls(true) {}

    virtual void spawn(Scene *scene, int num_start_samples, 
                       int max_bounce, 
                       BrdfPointLights &lights, bool is_photon = false) = 0;

    inline bool has_direct_vpls() const;
    inline void store_direct_vpls(bool direct);

protected:
    bool is_store_direct_vpls;
};

class VplGenerator : public IVplGenerator {
public:
    VplGenerator();

public:    
    /**
     * The generator assumes final gathering or VPL gathering to be performed, therefore, it traces 
     * and stores VPLs up to (max bounce - 1).
     *
     * Bounce is the number of path vertices excluding the vertex on the light and on the lens.
     */
    virtual void spawn(Scene *scene, int num_start_samples, int max_bounce, BrdfPointLights &lights, bool is_photon = false);

private:    
    virtual void generate_vpl(Scene *scene, int num_start_samples, int max_bounce, BrdfPointLights &vpls, bool is_photon = false);            
    void trace_vpl(Scene *scene, const Ray &r, Rgb T, int bounce, int max_bounce, BrdfPointLights &vpls, const BrdfPointLightRef &prev, Rgb path_L, Float path_pdf);    
}; 

inline bool IVplGenerator::has_direct_vpls() const {
    return is_store_direct_vpls;
}

inline void IVplGenerator::store_direct_vpls(bool direct) {
    this->is_store_direct_vpls = direct;
}

} // end namespace

#endif