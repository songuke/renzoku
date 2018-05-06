#ifndef _AGGREGATE_H_
#define _AGGREGATE_H_

#include "common.h"
#include "surface.h"

namespace Renzoku {
/**
 * Store a collection of surfaces.
 *
 * Inspired by Aggregate class in PBRT, but not a subclass of Surface
 * because Aggregate acts only as a wrapper to a collection of surfaces. 
 * It does not provide surface specific functions, e.g., BRDF sampling.
 */
class Aggregate {
public:
    Aggregate(const Surfaces &surfaces);
    virtual ~Aggregate();

    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    virtual bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;    

    /**
     * Visibility test between point P and Q. 
     */
    virtual bool hit(const Vec3 &p, const Vec3 &q, Float tmin, Float time) const;
    
    Float area() const;
    
    /**
     * Return a sample point from one of the surfaces in the scene. 
     * Uniform distribution.
     */
    void sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const;

    BoundingBox get_bounding_box() const;
    static BoundingBox get_bounding_box(Surfaces &surfaces);
    
    BoundingSphere get_bounding_sphere() const;

    /**
     * For visualization of Octree, BVH, etc.
     */
    virtual void get_bounding_boxes(BoundingBoxes &boxes) {}

    /**
     * Insert a new surface into spatial subdivision data structure (if any)
     */
    inline virtual void add_surface(Surface *surface);

    inline void set_scene(Scene *scene);

protected:
    /**
     * Compute probability density and cumulative density for surface sampling
     */
    void compute_density();
    
    /**
     * Compute and cache the scene AABB
     */
    void compute_bounding_box();
    void compute_bounding_sphere();
    

    /**
     * Set surface index to shape's data so that we can retrieve surface properties
     * after visibility test.
     */
    void set_surface_id();

    /**
     * Return the index of the shape we are going to sample
     * based on probability: A_j / sum(A_j). The probability of 
     * a point on a shape is 1 / A_j. So overall a point is 
     * sampled with probability 1 / sum(A_j), which is uniform.
     */    
    int sample_surface(Random &rd, Float &pdf) const;

protected:
    const Surfaces &surfaces;

    BoundingBox box;
    BoundingSphere sphere;

    Float *shape_pdf;
    Float *shape_cdf;
    Float shape_sum_area;

    Scene *scene;
}; 

inline void Aggregate::add_surface(Surface *surface) {
    // assume the just added surface is the last surface
    surface->get_shape()->set_surface_index(surfaces.size() - 1);
    compute_density();
}

inline void Aggregate::set_scene(Scene *scene) {
    this->scene = scene;
}

} // end namespace Renzoku

#endif

