#ifndef _BRDF_POINT_LIGHT_H_
#define _BRDF_POINT_LIGHT_H_

#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "rgb.h"
#include "light.h"
#include "material.h"
#include "onb.h"
#include "local_geometry.h"

namespace Renzoku {
     
enum BrdfPointLightTag {            // maintain for reproducing results in PG14.
    VPL_TAG_DEFAULT,
    VPL_TAG_NEW_GROUP,
    VPL_TAG_HEAL_ARTIFACTS,
    VPL_TAG_STANDARD_END,
    VPL_TAG_HEAL_END
};

struct BrdfPointLightRef {
    int index;

    BrdfPointLightRef(int index = -1) : index(index) {}

    bool is_valid() const { return index >= 0; }
};

class BrdfPointLight {
public:
    /**
     * Create a point light without material and parent point light. 
     * This indicates virtual point lights generated on the light sources.
     */ 
    virtual Light::Type get_light_type() const {
        return Light::BRDF_POINT_LIGHT;
    }

    BrdfPointLight() : distance(0.0f) {}
    
    /**
     * Create a VPL.
     *
     * The parent light can refer to a physical light, or another VPL (through index reference).
     */
    inline BrdfPointLight(const LocalGeometry &dg, Rgb contrib, Material *m, const Vec3 &wi, Light *ancestor_light);
    inline BrdfPointLight(const LocalGeometry &dg, Rgb contrib, Material *m, const Vec3 &wi, const BrdfPointLightRef &prev);    
    inline BrdfPointLight(const LocalGeometry &dg, Rgb power, Light *light);
        
    inline BrdfPointLightRef get_prev() const;
    inline Light *get_ancestor_light() const;
    inline void set_prev(const BrdfPointLightRef &ref);
    
    inline int get_bounce() const;
    inline void set_bounce(int bounce);

    inline Vec3 org() const;
    inline Vec3 normal() const;
    inline LocalGeometry get_dg() const;
    
    inline Rgb brdf(const Vec3 &wo) const;
    inline Rgb radiance(const Vec3 &wo) const;
	inline Rgb power() const;
        
	inline Float get_distance() const;
	inline void set_distance(Float distance);

    /**
     * Contribution = throughput / pdf. Contribution is equivalent to power.
     */
    inline Rgb contribution() const;
    inline void set_contribution(const Rgb &phi);

    /**
     * Return the incident ray.
     */
    inline Vec3 get_wi() const;

    /**
     * Return the local coordinate frame at the VPL.
     */
    inline Onb get_uvn() const;
    
    inline Material *get_material() const;
            
    inline Float get_radius() const;
    inline void  set_radius(Float radius);
    
    inline void set_shape(Shape *shape);
    inline Shape *get_shape() const;
    
    inline Vec2 get_texture_coords() const;

    /**
     * Throughput of the path. Does not include probability.
     */
    inline Rgb get_throughput() const;
    inline void set_throughput(const Rgb &T);

    inline void set_pdf(Float pdf);
    inline Float pdf() const;

    /**
     * Return true if VPL is on a light surface (represent direct lighting)
     */
    inline bool is_on_light() const;
        
    /**
     * Return true if VPL represents indirect lighting
     */ 
    inline bool is_indirect_light() const;

protected:    
    LocalGeometry dg;
    Vec3 wi;
    Rgb contrib;
    Material *m;    
    
    BrdfPointLightRef prev;     /// index of the previous VPL on the light subpath
    Light *ancestor_light;      /// mutual exclusive with prev

    int bounce;                 /// indicate the type of the VPL: bounce-0 means direct illumination.
    
    Float radius;               /// radius of the disk that estimates the area the VPL covers on the surface it is located on
        
    Shape *shape;               /// the shape on which the VPL is deposit. 
    
    Rgb throughput;
    Float path_pdf;

	Float distance;				/// distance to the previous path vertex

protected:
    BrdfPointLightTag tag;
 
public:
    inline BrdfPointLightTag get_tag() const;
    inline void set_tag(BrdfPointLightTag tag);
};

inline bool BrdfPointLight::is_on_light() const {
    return (m == NULL && ancestor_light != NULL);
}

inline bool BrdfPointLight::is_indirect_light() const {
    return (m != NULL && (prev.is_valid() || ancestor_light != NULL));
}

inline Vec3 BrdfPointLight::org() const { 
    return dg.p; 
}

inline Vec3 BrdfPointLight::normal() const { 
    return dg.n; 
}

inline LocalGeometry BrdfPointLight::get_dg() const {
    return dg;
}

inline Rgb BrdfPointLight::power() const {
    return contrib;
}

inline Rgb BrdfPointLight::contribution() const {
    return contrib;
}

inline void BrdfPointLight::set_contribution(const Rgb &phi) {
    contrib = phi;
}

inline Vec3 BrdfPointLight::get_wi() const {
    return wi;
}

inline Float BrdfPointLight::get_distance() const {
	return distance;
}

inline void BrdfPointLight::set_distance(Float distance) {
	this->distance = distance;
}

inline Onb BrdfPointLight::get_uvn() const {
    return dg.uvn;
}

inline Vec2 BrdfPointLight::get_texture_coords() const {
    return dg.uv;
}

inline BrdfPointLightRef BrdfPointLight::get_prev() const {    
    return prev;
}

inline void BrdfPointLight::set_prev(const BrdfPointLightRef &ref) {
    prev = ref;
}

inline Light *BrdfPointLight::get_ancestor_light() const {
    return ancestor_light;
}

inline int BrdfPointLight::get_bounce() const {
    return bounce;
}

inline void BrdfPointLight::set_bounce(int bounce) {
    this->bounce = bounce;
}

inline Material *BrdfPointLight::get_material() const {
    return m;
}

inline Float BrdfPointLight::get_radius() const {
    return radius;
}

inline void  BrdfPointLight::set_radius(Float radius) {
    this->radius = radius;
}

inline Shape *BrdfPointLight::get_shape() const {
    return shape;
}

inline void BrdfPointLight::set_shape(Shape *shape) {
    this->shape = shape;
}

inline Rgb BrdfPointLight::brdf(const Vec3 &wo) const {
    if (this->is_on_light())
        return DefaultRgb::white;

    return m->eval(dg, wi, wo);    
}

inline void BrdfPointLight::set_pdf(Float pdf) {
    this->path_pdf = pdf;
}

inline Float BrdfPointLight::pdf() const {
    return path_pdf;
}

inline void BrdfPointLight::set_throughput(const Rgb &T) {
    this->throughput = T;
}

inline Rgb BrdfPointLight::get_throughput() const {
    return throughput;
}

inline BrdfPointLightTag BrdfPointLight::get_tag() const {
    return tag;
}
    
inline void BrdfPointLight::set_tag(BrdfPointLightTag tag) {
    this->tag = tag;
}

inline BrdfPointLight::BrdfPointLight(const LocalGeometry &dg, Rgb power, Material *m, const Vec3 &wi, const BrdfPointLightRef &prev) 
    : dg(dg), contrib(power), m(m), prev(prev), wi(wi), ancestor_light(NULL), distance(0.0f)
{    
    
}

inline BrdfPointLight::BrdfPointLight(const LocalGeometry &dg, Rgb power, Material *m, const Vec3 &wi, Light *ancestor_light) 
    : dg(dg), contrib(power), m(m), prev(BrdfPointLightRef()), wi(wi), ancestor_light(ancestor_light), distance(0.0f)
{
    
} 

inline BrdfPointLight::BrdfPointLight(const LocalGeometry &dg, Rgb power, Light *light) 
    : dg(dg), contrib(power), m(NULL), prev(BrdfPointLightRef()), wi(Vec3()), ancestor_light(ancestor_light), distance(0.0f)
{
    
} 

inline Rgb BrdfPointLight::radiance(const Vec3 &wo) const {    
    Rgb brdf = this->brdf(wo);        
	return contrib * brdf;
}

} // end namespace

#endif

