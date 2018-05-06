#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "shape.h"
#include "material.h"
#include "gl_vertex.h"
#include "sphere.h"

#include <vector>
using namespace std;

namespace Renzoku {

class AreaLight;

/**
 * Represent a surface that consists of a shape and a material. A surface can also be an area light source.
 *
 * Surface is also responsible to flipping shading normals, and establishing local frame at the intersection of a ray-shape pair. * 
 * The shading normal can be flipped in order to make sure it points in the same direction with wo.
 *
 * The stored normals of shapes are not flipped, and can be queried independently.
 */ 
class Surface {
public:
    /**
     * Create a normal surface with a material
     */
    Surface(Shape *_shape, Material *_material);

    /**
     * Create a surface and attach an area light to it
     */
    Surface(Shape *_shape, const Rgb &emission);
    Surface(Shape *_shape, AreaLight *light);

    /**
     * Create a sphere and attach an environment light to it
     */
    Surface(Sphere *_sphere, ImageFloat *envmap);

    Surface(const Surface &s);
        
    bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;
        
    /**
     * Geometry/shape sampling interface. Depending on which type of property 
     * is attached to a surface, such sampling might or might not be used.
     */
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;         
    virtual Float pdf(const Vec3 &p) const;

    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const;         
    virtual Float pdf(const Vec3 &p, const Receiver &patch) const;

    inline void set_id(ID id);
    inline virtual ID get_id() const;

    inline Float        area() const;
    inline Material*    get_material() const;
    inline AreaLight*   get_area_light() const;
    inline EnvLight*    get_env_light() const;
    inline Shape*       get_shape() const;

    inline BoundingBox      get_bounding_box() const;
    inline BoundingSphere   get_bounding_sphere() const;

    /**
     * Fill the data structure with all vertex coordinates, normals, and materials. For OpenGL rendering. 
     */
    virtual void get_vertex_data(vector<GlVertex> &vertices); 

    /**
     * Return true if this is an area light.
     */
    inline bool is_light() const;

    /** 
     * Return true if this is an environment light.
     */
    inline bool is_env_light() const;

    inline void reverse_normal(bool flipped);
    inline bool is_reverse_normal() const;

    /** 
     * Return a local orthogonal basis
     */
    inline Onb get_basis(const Vec3 &p) const;

    inline bool get_two_sided() const;
    inline void set_two_sided(bool two);

    void fill_hit_record(const Ray &r, HitRecord &hit) const;

protected: 
    ID id;
    Shape *shape;
    
    Material *material;
    AreaLight *light;
    EnvLight *env_light;

    // Global basis can only be maintained for flat shape. Not suitable for parametric shape. 
    // Also, basis depends on flip status. 
    //Onb uvn;                // global orthogonal basis for the shape, for consistent wi approximation
    
    bool is_flipped;
    bool is_two_sided;    
};

class Surfaces : public vector<Surface> {

private:  // force to use append
    void push_back(const Surface &s) {
        vector<Surface>::push_back(s);
    }
    void push_back(Surface &s) {
        vector<Surface>::push_back(s);
    }

public:
    void append(const Surface &s) {
        s.get_shape()->set_surface_index(this->size());
        this->push_back(s);
    }
};

inline Float Surface::area() const {
    return shape->area();
}

inline Material* Surface::get_material() const { 
    return material; 
}

inline AreaLight* Surface::get_area_light() const {
    return light;
}

inline EnvLight *Surface::get_env_light() const {
    return env_light;
}

inline bool Surface::is_light() const {
    return light != NULL;
}

inline bool Surface::is_env_light() const {
    return env_light != NULL;
}

inline Shape* Surface::get_shape() const { 
    return shape; 
}

inline BoundingBox Surface::get_bounding_box() const { 
    return shape->get_bounding_box(); 
}

inline BoundingSphere Surface::get_bounding_sphere() const {
    return shape->get_bounding_sphere();
}

inline void Surface::reverse_normal(bool flipped) {
    is_flipped = flipped;
}

inline Onb Surface::get_basis(const Vec3 &p) const {
    Vec3 n = shape->normal(p);
    Vec3 u = shape->tangent();
    Onb uvn;
    uvn.init_from_nu(n, u);
    return uvn;
}

inline bool Surface::is_reverse_normal() const {
    return is_flipped;
}

inline void Surface::set_id(ID id) {
    this->id = id;
}

inline ID Surface::get_id() const {
    return id;
}

inline void Surface::set_two_sided(bool two) {
    is_two_sided = two;
}

inline bool Surface::get_two_sided() const {
    return is_two_sided;
}

} // end namespace Renzoku


#endif

