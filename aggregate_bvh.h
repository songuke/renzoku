#ifndef _AGGREGATE_BVH_H_
#define _AGGREGATE_BVH_H_

#include "common.h"
#include "aggregate.h"
#include "boundingbox.h"

namespace Renzoku {

struct BvhNode {
    BoundingBox box;
    int start, end;                 // index to the sorted morton array
    int index_left, index_right;    
    int level;
};

/** 
 * CPU implementation of linear BVH.
 */
class AggregateBvh : public Aggregate {
public:
    AggregateBvh(Surfaces &surfaces);
    ~AggregateBvh();

    inline bool hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const;
    inline bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;

	inline void get_bounding_boxes(BoundingBoxes &boxes);
    
    inline virtual void add_surface(Surface *surface);
        
protected:    
    void build();
    void cleanup();

protected:        
    BvhNode *root;    
    BvhNode *nodes;
    int num_nodes;
    int num_levels;
        
    Triangle            *sorted_triangles;        
    unsigned int        *sorted_morton;    
    int                 *sorted_surface_id;
    BoundingBox         *sorted_boxes;
    int                 num_triangles;    
}; 

inline void AggregateBvh::get_bounding_boxes(BoundingBoxes &boxes) {
	boxes.clear();
    for (int i = 0; i < num_nodes; ++i) {
        if (nodes[i].level == 2)
            boxes.push_back(nodes[i].box);
    }
}

inline void AggregateBvh::add_surface(Surface *surface) {
    Aggregate::add_surface(surface);
    // not support incremental build yet
    // simply rebuild all
    cleanup();
    build();    
}

} // end namespace Renzoku

#endif

