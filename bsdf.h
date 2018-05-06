#ifndef _BSDF_H_
#define _BSDF_H_

namespace Renzoku {

class Bsdf {
public:
    enum Type {
        LAMBERTIAN = 0,
        WARD = 1,
        PHONG = 2,
        GLASS = 3,
        MIRROR = 4,
        THIN_TRANSPARENT = 5,
        THIN_TRANSLUCENT = 6
    };
    
    virtual Bsdf::Type get_bsdf_type() const = 0;

    /**
     * Sample an in-going direction from the BRDF given an out-going direction.
     */
    virtual Rgb sample(Random &rd, const Onb &uvn, const Vec3 &wo, Vec3 &wi, Float &pdf) = 0;
        
    /**
     * Evaluate the BRDF value given an in, out direction. 
     */
    virtual Rgb eval(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) = 0;
    
    /**
     * Compute the probability of sampling the direction \wi.
     */
    virtual Float pdf(const Onb &uvn, const Vec3 &wo, const Vec3 &wi) = 0;
    
    /** 
     * Return true if there is only a single ray that can be sampled from this material. Example: glass, mirror.
     */
    virtual bool is_singular() const { return false; }

    /**
     * Return true if the material lets ray pass through.
     */
    virtual bool is_transmissive() const { return false; }

    bool is_specular() const {
        return this->get_bsdf_type() != LAMBERTIAN;
    }

    bool is_diffuse() const {
        return this->get_bsdf_type() == LAMBERTIAN;
    }

};

class IDiffuse { 
public:
    virtual Rgb get_diffuse_component() const = 0;
};

class ISpecular {
public:
    virtual Rgb get_specular_component() const = 0;

    /**
     * The degree of glossiness approximated by the Phong exponent.
     *
     * This is for classifying glossy materials, not for sampling or evaluation.
     */
    virtual Float get_gloss_exponential() const = 0;
};

} // end namespace

#endif