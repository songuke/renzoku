#ifndef _SHAPE_H_
#define _SHAPE_H_

#include "common.h"
#include "named_object.h"
#include "base_object.h"
#include "vec2.h"
#include "vec3.h"
#include "ray.h"
#include "onb.h"
#include "random.h"
#include "boundingbox.h"
#include "gl_vertex.h"
#include "hitrecord.h"

namespace Renzoku {

/**
 * \brief Shape represents basic geometry objects. It only contains vertex coordinates. 
 * Shape can be extended to be a set of shapes which groups similar shapes with same materials together under a surface.
 */
class Shape : public BaseObject, public NamedObject {
public:
    Shape();    

    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const = 0;
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const = 0;

    virtual Float area() const = 0;
    
    /**
     * Sample a point on the shape surface given the random parameters.
     */
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const = 0; 
    
    /**
     * Probability of a point being generated on the surface of the shape.
     */
    virtual Float pdf(const Vec3 &p) const = 0;

    /**
     * Sample a point on the shape surface depending on properties of the receiving patch such as location and orientation.
     * 
     * Useful for sampling points on light sources. 
     */
    virtual void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf, const Receiver &patch) const {
        this->sample(rd, p, n, pdf);
    }

    /**
     * Probability of the point generated with receiving patch considered.
     */
    virtual Float pdf(const Vec3 &p, const Receiver &patch) const {
        return this->pdf(p);
    }

    /**
     * Add vertex coordinates to the array. 
     */        
    virtual void get_triangle_coordinates(vector<Float> &coords) const = 0;
    virtual Triangle *get_triangle_pointer(int triangle_index) const { return NULL; }

    /**
     * Add vertex coordinates and normals to the array. 
     */ 
    virtual void get_vertex_data(vector<GlVertex> &vertices) const = 0;

    virtual void get_vertex_positions(vector<Vec3> &positions) const = 0;
    virtual void set_vertex_normal(int vertex_idx, const Vec3 &normal) = 0;

    /**
     * Return the number of triangles that this shape represents.
     */
    virtual int get_triangle_count() const { return 0; }
    virtual void get_triangles(vector<Triangle> &list) const {}

    /**
     * Return the axis-aligned bounding box of the shape
     */
	BoundingBox get_bounding_box() const {
        return box;
    }

    BoundingSphere get_bounding_sphere() const {
        return sphere;
    }

    /**
     * Facet normal
     */
    virtual Vec3 normal(const Vec3 &p) const = 0;

    /**
     * Facet tangent
     */
    virtual Vec3 tangent() const = 0;

    /**
     * Shading normals can be interpolated from vertex normals, 
     * or obtained from bump maps.
     */
    virtual Vec3 shading_normal(const Vec3 &p) const = 0;

    /**
     * Return interpolated UV coordinates.
     */ 
    virtual Vec2 texture_uv(const Vec3 &p) const = 0;

    inline int get_surface_index() const;
    inline virtual void set_surface_index(int index);

    inline int get_smooth_group() const;
    inline void set_smooth_group(int idx);
    
    inline int get_object_index() const;
    inline void set_object_index(int idx);

    /**
     * Fill remaining parameters of hit record based on t, beta, and gamma.
     */
    virtual void fill_hit_record(const Ray &r, const GeometryHit &gh, HitRecord &record) const {
        record.copy_geometry_hit(gh);        
    }
    
protected:
    virtual void compute_bounding_box() = 0;
    virtual void compute_bounding_sphere() = 0;

    BoundingBox box;
    BoundingSphere sphere;

protected:
    int surface_index;         // for looking up in the surface list

    int smooth_group;    
    int object_index;          // for controlling material of geometric elements of an object
};

inline Shape::Shape() : surface_index(-1), smooth_group(-1), object_index(-1) {
    
}

inline int Shape::get_smooth_group() const {
    return smooth_group;
}

void Shape::set_smooth_group(int idx) {
    smooth_group = idx;
}

inline int Shape::get_surface_index() const {
    return surface_index;
}

inline void Shape::set_surface_index(int idx) {
    surface_index = idx;
}

inline void Shape::set_object_index(int idx) {
    object_index = idx;
}

inline int Shape::get_object_index() const {
    return object_index;
}

};
#endif

