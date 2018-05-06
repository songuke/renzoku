#include "surface.h"
#include "log.h"
#include "ward.h"

#include "area_light.h"
#include "env_light.h"
#include "glass.h"

#include "resource.h"

namespace Renzoku {

Surface::Surface(Shape *_shape, Material *_material) 
    : shape(_shape), material(_material), light(NULL), env_light(NULL), is_flipped(false), is_two_sided(false) {
    
    set_id(Resource::get_next_id());
}

Surface::Surface(Shape *_shape, const Rgb &emission)
    : shape(_shape), material(NULL), env_light(NULL), is_flipped(false), is_two_sided(false) {

    set_id(Resource::get_next_id());

    light = new AreaLight(this, emission);
}

Surface::Surface(Shape *_shape, AreaLight *light) : shape(_shape), light(light), is_flipped(false), is_two_sided(false) {
    set_id(Resource::get_next_id());
}

Surface::Surface(const Surface &s)
    : shape(s.shape), material(s.material), light(s.light), env_light(s.env_light),
      is_flipped(s.is_flipped), is_two_sided(s.is_two_sided)
{
    if (light)
        light->set_surface(this);
}

Surface::Surface(Sphere *_sphere, ImageFloat *envmap) 
    : shape(_sphere), material(NULL), light(NULL), is_flipped(true) /* for inner sphere */, is_two_sided(false) {
    
    set_id(Resource::get_next_id());

    env_light = new EnvLight(this, envmap);
}

bool Surface::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    GeometryHit gh;
    shape->hit(r, tmin, tmax, time, gh);
    record.copy_geometry_hit(gh);        
    if (gh.hit) {
        record.shape->fill_hit_record(r, gh, record);        
        this->fill_hit_record(r, record);
    }
    return gh.hit;
}

bool Surface::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return shape->hit(r, tmin, tmax, time);
}

void Surface::fill_hit_record(const Ray &r, HitRecord &record) const {
    record.material = material;
    record.light = light;
    record.env_light = env_light;
    record.shape = shape;
    record.surface = (Surface *)this;

    // displacement map can overwrite shading normals interpolated from vertex normals
    /*
    if (material && material->has_displacement_map()) {
        Texture *disp = material->get_displacement_map();
        Vec3 d = disp->lookup(record.uv).to_vec3();
        d = unit_vector(Vec3(d.x() * 2.0f - 1.0f, d.y() * 2.0f - 1.0f, d.z()));
        record.shading_normal = record.uvn->local_to_world(d);
    } else {
        record.shading_normal = record.normal;
    }*/

    // NOTE: this flipped is not permanent and is not reflected in OpenGL get_vertex_data.
    if (is_flipped) { // due to vertex ordering
        record.normal = -record.normal;
        record.shading_normal = -record.shading_normal;
    }
    
    if (is_two_sided && dot(record.shading_normal, r.dir()) > 0.0f) {
        if (record.material && record.material->get_bsdf()->get_bsdf_type() != Bsdf::GLASS) {  // glass cannot be two sided
            record.shading_normal = -record.shading_normal; 
            record.normal = -record.normal;
        }
    }

    // make basis after handling normal flipping
    // this ensures correctly flipped basis is used for BRDF sampling in integrators.
    record.uvn.init_from_nu(record.shading_normal, record.tangent);
}
    
void Surface::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    shape->sample(rd, p, n, pdf);
    if (is_flipped)
        n = -n;
}

Float Surface::pdf(const Vec3 &p) const {
    return shape->pdf(p);
}

void Surface::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const {
    shape->sample(rd, p, n, pdf, patch);
    if (is_flipped)
        n = -n;
}

Float Surface::pdf(const Vec3 &p, const Receiver &patch) const {
    return shape->pdf(p, patch);
}

void Surface::get_vertex_data(vector<GlVertex> &vertices) {
    // for environment light, skip
    // for area light, assume it is black since area light can block light from other surfaces.
    if (env_light != NULL) return;      
    
    int old_size = vertices.size();
    shape->get_vertex_data(vertices);

    for (int i = old_size; i < vertices.size(); ++i) {
        GlVertex *v = &vertices[i];
                
        if (light) {
            v->material_type = Bsdf::LAMBERTIAN;
                        
            v->kd[0] = 0.0f; 
            v->kd[1] = 0.0f; 
            v->kd[2] = 0.0f; 
                        
            v->ks[0] = 0.0f; 
            v->ks[1] = 0.0f; 
            v->ks[2] = 0.0f; 
                        
            v->glossy[0] = 1.0f;
            v->glossy[1] = 1.0f;

        } else {
            
            MarshalObject mo;
            material->marshal(mo);                        
            
            v->material_type = material->get_bsdf()->get_bsdf_type();

            Rgb kd = mo.dict.get_vec3("kd", Vec3(0.0f));
            v->kd[0] = kd.red(); 
            v->kd[1] = kd.green();
            v->kd[2] = kd.blue();

            Rgb ks = mo.dict.get_vec3("ks", Vec3(0.0f));
            v->ks[0] = ks.red();
            v->ks[1] = ks.green();
            v->ks[2] = ks.blue();

            Vec2 glossy = mo.dict.get_vec2("glossy", Vec2(0.0f));
            v->glossy[0] = glossy.x();
            v->glossy[1] = glossy.y();            
        }
    }
}

} // end namespace Renzoku
