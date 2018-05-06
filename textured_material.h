#ifndef _TEXTURED_MATERIAL_H_
#define _TEXTURED_MATERIAL_H_

#include "material.h"
#include "texture.h"

namespace Renzoku {

/**
 * A decorator that adds texture lookup to BRDF.
 */
class TexturedMaterial : public Material {
public:
    TexturedMaterial(Material *m, Texture *t) : m(m), tex(t), tex_disp(NULL) {
    }

    virtual Material::Type get_material_type() const {
        return Material::TEXTURED;
    }
    
    virtual void marshal(MarshalObject &mo) {
        m->marshal(mo);
    }

    bool has_texture() const { return tex != NULL; }
    bool has_displacement_map() const { return tex_disp != NULL; }

    Texture *get_texture() const { return tex; }
    void set_texture(Texture *tex) { this->tex = tex; }

    Texture *get_displacement_map() const { return tex_disp; }
    void set_displacement_map(Texture *tex) { this->tex_disp = tex; }
        
    virtual Bsdf* get_bsdf() {
        return m->get_bsdf();
    }

    virtual Rgb sample(Random &rd, const LocalGeometry &dg, const Vec3 &wo, Vec3 &wi, Float &pdf) {
        Rgb brdf = m->sample(rd, dg, wo, wi, pdf);        
        if (tex && pdf > 0.0f && brdf != DefaultRgb::black) {
            return brdf * tex->lookup(dg.uv);
        } else return brdf;
    }
  
    virtual Rgb eval(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) {
        Rgb brdf = m->eval(dg, wo, wi);
        if (tex && brdf != DefaultRgb::black) {
            return brdf * tex->lookup(dg.uv);
        } else return brdf;
    }

    virtual Float pdf(const LocalGeometry &dg, const Vec3 &wo, const Vec3 &wi) {
        return m->pdf(dg, wo, wi);
    }

protected:
    Material *m;
    Texture *tex;
    Texture *tex_disp;
};

} // end namespace

#endif