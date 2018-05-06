#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "common.h"
#include "named_object.h"
#include "rgb.h"
#include "vec2.h"
#include "local_geometry.h"
#include "hitrecord.h"
#include "bsdf.h"
#include "marshal.h"

#include "lambertian.h"
#include "phong.h"
#include "ward.h"
#include "mirror.h"
#include "glass.h"

namespace Renzoku {

/**
 * Material class defines the common interface for all materials and material decorators. 
 *
 * The Material and BsdfMaterial have the same interface. However, we do not want to merge them
 * since TexturedMaterial (or other decorator materials) inherits from Material but it should not store a BSDF internally. 
 * It should rely on the material it decorates for BSDF sampling.
 */
class Material : public NamedObject, public IMarshalable {    
public:
    Material() : index(-1) {}

    /**
     * Fill in the general material that can be transferred to GPU or network.
     */
    virtual void marshal(MarshalObject &m) = 0;

    /**
     * Returns a color to represent this material. 
     * This is used in displaying the scene in wireframe or preview mode (non-physically based rendering).
     * 
     * White is returned by default.
     */
    virtual Rgb get_representative_color() const { return DefaultRgb::white; }

    virtual Bsdf* get_bsdf() = 0;
    virtual Rgb sample(Random &rd, const LocalGeometry &dg, const Vec3 &wo, Vec3 &wi, Float &pdf) = 0;
    virtual Rgb eval(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) = 0;
    virtual Float pdf(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) = 0;

    inline int     get_material_index() const;
    inline void    set_material_index(int index);
    
    enum Type {
        BSDF,
        TEXTURED
    };
    virtual Material::Type get_material_type() const = 0;

    bool has_texture() const {
        return this->get_material_type() == TEXTURED;
    }

protected:
    int index;              // index to the global material array
};

class BsdfMaterial : public Material {
public:
    BsdfMaterial(Bsdf *bsdf) : bsdf(bsdf) {}
    
    virtual Material::Type get_material_type() const {
        return Material::BSDF;
    }

    virtual void marshal(MarshalObject &m) {
        m.dict.clear();

        // for consistency, we use a switch here. Another option is to implement marshal in each BSDF,
        // but we don't want to pollute such classes.
        m.dict.add("bsdf_type", (int)bsdf->get_bsdf_type());

        switch (bsdf->get_bsdf_type()) {
            case Bsdf::LAMBERTIAN:
            {
                m.dict.add("kd", ((Lambertian *)bsdf)->get_diffuse_component().to_vec3());
                break;
            }
            case Bsdf::PHONG: 
            {
                ModifiedPhong *phong = (ModifiedPhong *)bsdf;
                m.dict.add("kd", phong->get_diffuse_component().to_vec3());
                m.dict.add("ks", phong->get_specular_component().to_vec3());
                m.dict.add("glossy", phong->get_glossy_coefficients());
                break;                          
            }
            case Bsdf::WARD:
            {
                Ward *ward = (Ward *)bsdf;
                m.dict.add("kd", ward->get_diffuse_component().to_vec3());
                m.dict.add("ks", ward->get_specular_component().to_vec3());
                m.dict.add("glossy", ward->get_glossy_coefficients());
                break;
            }
        }
    }

    /**
     * Sample an in-going direction from the BRDF given an out-going direction.
     */
    virtual Rgb sample(Random &rd, const LocalGeometry &dg, const Vec3 &wo, Vec3 &wi, Float &pdf) {
        return bsdf->sample(rd, dg.uvn, wo, wi, pdf);
    }
        
    /**
     * Evaluate the BRDF value given an in, out direction. 
     */
    virtual Rgb eval(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) {
        return bsdf->eval(dg.uvn, wo, wi);
    }
    
    /**
     * Compute the probability of sampling the direction \wi.
     */
    virtual Float pdf(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) {
        return bsdf->pdf(dg.uvn, wo, wi);
    }
    
    /**
     * Return a unique ID for the material type
     */
    virtual Bsdf* get_bsdf() { return bsdf; }
    
protected:    
    Bsdf *bsdf;    
};

int Material::get_material_index() const {
    return index;
}

void Material::set_material_index(int index) {
    this->index = index;
}

class DefaultMaterial {
public:
    static BsdfMaterial *white();
    static BsdfMaterial *pi();

protected:
    static BsdfMaterial *_white;
    static BsdfMaterial *_pi;
};

} // end namespace

#endif

