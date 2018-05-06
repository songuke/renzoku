#ifndef _OCTREE_H_
#define _OCTREE_H_

#include "vec3.h"
#include "boundingbox.h"
#include "surface.h"

#include <cassert>
#include <algorithm>
using namespace std;

namespace Renzoku{

struct OctreeNode {
	BoundingBox box;            /// the node's bounding box    
    Shapes shapes;               /// the list of shapes that OctreeNode contains; could be empty.
    OctreeNode *child[8];

	OctreeNode(const BoundingBox &bb) : box(bb) {		
        for (int i = 0; i < 8; ++i) child[i] = NULL;
	}

    ~OctreeNode() {
        if (! is_leaf())
            for (int i = 0; i < 8; ++i) delete child[i]; // call destructors of child recursively        
    }

    /**
     * Leaf node has no children.
     */
    bool is_leaf() const {
        /*
        for (int i = 0; i < 8; ++i) {
            if (child[i]) return false;
        }
        return true;*/

        // since we either have 8 child nodes or none, just need to check child[0]
        return child[0] == NULL;
    }

    /**
     * Split the bounding box into 8 smaller boxes and attach to 8 child nodes.
     * A proper child is selected and the current surface is inserted to that child node.
     */
    void split() {
        Vec3 m = box.v_min;
        Vec3 M = box.v_max;
        Vec3 centroid = (m + M) * 0.5;
        Vec3 axis = centroid - m;

        child[0] = new OctreeNode(BoundingBox(m,                               centroid));
        child[1] = new OctreeNode(BoundingBox(m + Vec3(0, 0,        axis.z()), centroid + Vec3(0, 0,        axis.z())));
        child[2] = new OctreeNode(BoundingBox(m + Vec3(0, axis.y(), 0),        centroid + Vec3(0, axis.y(), 0)));
        child[3] = new OctreeNode(BoundingBox(m + Vec3(0, axis.y(), axis.z()), centroid + Vec3(0, axis.y(), axis.z())));

        child[4] = new OctreeNode(BoundingBox(m + Vec3(axis.x(), 0, 0),               centroid + Vec3(axis.x(), 0, 0)));
        child[5] = new OctreeNode(BoundingBox(m + Vec3(axis.x(), 0,        axis.z()), centroid + Vec3(axis.x(), 0,        axis.z())));
        child[6] = new OctreeNode(BoundingBox(m + Vec3(axis.x(), axis.y(), 0),        centroid + Vec3(axis.x(), axis.y(), 0)));
        child[7] = new OctreeNode(BoundingBox(m + Vec3(axis.x(), axis.y(), axis.z()), centroid + Vec3(axis.x(), axis.y(), axis.z())));

        /*
        child[4] = new OctreeNode(new BoundingBox(centroid + Vec3(0, -axis.y(), -axis.z()), M + Vec3(0, -axis.y(), -axis.z())));
        child[5] = new OctreeNode(new BoundingBox(centroid + Vec3(0, -axis.y(),  0),        M + Vec3(0, -axis.y(),  0)));
        child[6] = new OctreeNode(new BoundingBox(centroid + Vec3(0,  0,        -axis.z()), M + Vec3(0,  0,        -axis.z())));
        child[7] = new OctreeNode(new BoundingBox(centroid,                                 M));
        */
        
        // re-insert
        //assert(shapes.size() == 1); // guarantee that when split there is only a surface in the node
        for (int k = 0; k < shapes.size(); ++k) {
            BoundingBox shape_box = shapes[k]->get_bounding_box();
            for (int i = 0; i < 8; ++i) {
                if (child[i]->box.overlap(&shape_box))
                    child[i]->shapes.push_back(shapes[k]);
            }
        }
        this->shapes.clear();
    }

    /**
     * When a node that has all child nodes having similar set of surfaces, that means previous split is useless. 
     * Take the surface set and attach to the parent and delete child nodes. 
     */
    void condense() {
        if (is_leaf()) return;

        // all child must be leaf nodes to condense
        for (int i = 0; i < 8; ++i) 
            if (! child[i]->is_leaf()) return;

        for (int i = 1; i < 8; ++i) {            
            if (child[i]->shapes.size() != child[i - 1]->shapes.size())
                return;

            for (int k = 0; k < child[i]->shapes.size(); ++k) {                
                if (child[i]->shapes[k] != child[i - 1]->shapes[k])
                    return;
            }
        } 

        // transfer the set into parent node
        shapes = child[0]->shapes;
                
        for (int i = 0; i < 8; ++i) {
            delete child[i];
            child[i] = NULL;
        }        
    }
    
    /**
     * Sort the node shapes according to pointer address.
     */
    void sort() {        
        std::sort(shapes.begin(), shapes.end());
    }
};

/**
 * A simple structure for collecting Octree statistics.
 */
struct OctreeStats {
    int num_surfaces;
    int num_nodes;
    int num_leaves;

    int thinnest_node;                      // only consider leaf nodes
    int thinnest_node_num_surfaces;
    
    int fattest_node;
    int fattest_node_num_surfaces;

    OctreeStats() {
        num_surfaces = 0;
        num_nodes = 0;
        num_leaves = 0;

        thinnest_node = -1;        
        thinnest_node_num_surfaces = 1e9;

        fattest_node = -1;
        fattest_node_num_surfaces = 0;
    }
};

class Octree {
protected:
    OctreeNode *root;
    int leaf_capacity;      /// the number of surfaces a leaf node can store before splitting occurs. Default: 4.
    int max_level;          /// maximum level (height) of the tree. Default: 8.

public:
    Octree(const BoundingBox &box);

    ~Octree() {
        if (root) delete root;
    }
    
    void insert(Shape *s);
    
    void condense();

    /**
     * Check correctness of the octree: non-leaf node should not contain any surfaces.
     * Also collect internal octree stats.
     */
    void verify();

    /**
     * Define how many surfaces a leaf node can store.
     *
     * The smaller this value is, the more tree nodes are generated and can cause overhead during ray-octree intersection.
     *
     * Default: 16.
     */
    void set_leaf_capacity(int size);

    /**
     * Define the maximum height of the tree. When the height of the tree is reached during insertion, all surfaces are
     * simply inserted into the leaf node. 
     *
     * Default: 8.
     */
    void set_max_level(int num);

	/**
	 * Return all bounding boxes in the tree. For debugging purpose.
	 */
	void get_bounding_boxes(BoundingBoxes &boxes);

    bool hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const;
    bool hit(const Ray &r, Float tmin, Float tmax, Float time) const;
 
    
protected:
    void insert_node(OctreeNode *node, Shape *s, int level, int max_level);

    void condense_node(OctreeNode *node);

	void get_bounding_boxes_node(OctreeNode *node, BoundingBoxes &boxes);

    bool hit_node(OctreeNode *node, const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const;
    bool hit_node(OctreeNode *node, const Ray &r, Float tmin, Float tmax, Float time) const;

    void verify_node(OctreeNode *node, OctreeStats &stats); 
};

};
#endif
