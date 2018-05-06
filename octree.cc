#include "octree.h"
#include "log.h"
#include <cassert>

namespace Renzoku {

Octree::Octree(const BoundingBox &box) {
    root = new OctreeNode(box);
    leaf_capacity = 16;
    max_level = 8;
}

void Octree::set_leaf_capacity(int size) {
    leaf_capacity = size;
}

void Octree::set_max_level(int num) {
    max_level = num;
}

void Octree::insert_node(OctreeNode *node, Shape *s, int level, int max_level) {    
    if (level == max_level) {
        node->shapes.push_back(s);
        return;
    }
        
    if (node->is_leaf()) {
        if (node->shapes.size() < leaf_capacity) {
            node->shapes.push_back(s);
            return;
        } else { // leaf node occupied then split (spawn 8 child nodes)
            node->split();
        }
    }
    
    // find the child nodes that surface can be inserted. 
    // If the surface spreads over bounding boxes of more than one child node, 
    // it will be inserted to all of the child nodes that the surface overlaps.
    BoundingBox box = s->get_bounding_box();
    for (int i = 0; i < 8; ++i) {
        if (node->child[i]->box.overlap(&box)) {
            insert_node(node->child[i], s, level + 1, max_level);
        }
    }
}

void Octree::insert(Shape *s) {
    insert_node(root, s, 1, max_level); // (8^8 - 1) / (8 - 1) = 2^24 / 7 ~ 2 million nodes in total
}

void Octree::condense_node(OctreeNode *node) {
    if (node->is_leaf()) {        
        node->sort();
        return;
    }

    for (int i = 0; i < 8; ++i) {
        condense_node(node->child[i]);
    }
    node->condense();
}

void Octree::condense() {
    condense_node(root);
}

bool Octree::hit(const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const {
    return hit_node(root, r, tmin, tmax, time, record);    
}

bool Octree::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    return hit_node(root, r, tmin, tmax, time);
}

bool Octree::hit_node(OctreeNode *node, const Ray &r, Float tmin, Float tmax, Float time, GeometryHit &record) const {    
    if (! node->box.hit(r, tmin, tmax)) {
        //cout << "Ray miss the box." << endl;
        //exit(0);
        return false;
    }
    if (node->is_leaf()) {
        // brute-force
        bool hit = false;
        for (int k = 0; k < node->shapes.size(); ++k) {
            GeometryHit gh;
            if (node->shapes[k]->hit(r, tmin, tmax, time, gh)) {
                tmax = gh.t;
                record = gh;
                hit = true;                
            }
        }
        return hit;
    } else {
        // this ray can hit more than one child bounding boxes, check all boxes. Worst case
        // this ray hits the common point of all eight bounding boxes (the centroid).         
        bool hit = false;
        for (int i = 0; i < 8; ++i) {
            GeometryHit gh;
            if (hit_node(node->child[i], r, tmin, tmax, time, gh)) {
                tmax = gh.t;
                record = gh;
                hit = true;
            }
        }
        return hit;
    }
}

bool Octree::hit_node(OctreeNode *node, const Ray &r, Float tmin, Float tmax, Float time) const {
    if (! node->box.hit(r, tmin, tmax)) return false;

    if (node->is_leaf()) {
        // brute-force
        for (int k = 0; k < node->shapes.size(); ++k) {
            if (node->shapes[k]->hit(r, tmin, tmax, time)) {
                return true;
            }
        }
        return false;
    } else {
        for (int i = 0; i < 8; ++i) {
            if (hit_node(node->child[i], r, tmin, tmax, time)) {
                return true;
            }
        }
        return false;
    }
}

void Octree::verify_node(OctreeNode *node, OctreeStats &stats) {
    if (node == NULL) return;
    ++stats.num_nodes;

    if (! node->is_leaf()) {
        if (node->shapes.size() > 0)
            Log::error() << "Non-leaf contains shapes!" << endn;
    } else {
        ++stats.num_leaves;

        if (node->shapes.size() > 0) {
            stats.num_surfaces += node->shapes.size();
            //cout << "Found leaf with shapes: " << node->shapes.size() << endl;
        }

        if ((int)node->shapes.size() < stats.thinnest_node_num_surfaces)
            stats.thinnest_node_num_surfaces = node->shapes.size();

        if ((int)node->shapes.size() > stats.fattest_node_num_surfaces)
            stats.fattest_node_num_surfaces = node->shapes.size();
    }

    for (int i = 0; i < 8; ++i)
        verify_node(node->child[i], stats);
}

void Octree::verify() {
    OctreeStats stats;
    verify_node(root, stats);

    Log::info() << "Octree surfaces        : " << stats.num_surfaces << endn;
    Log::info() << "Octree nodes           : " << stats.num_nodes << endn;
    Log::info() << "Octree leaves          : " << stats.num_leaves << endn;
    Log::info() << "Thinnest node surfaces : " << stats.thinnest_node_num_surfaces << endn;
    Log::info() << "Fattest node surfaces  : " << stats.fattest_node_num_surfaces << endn;    
}

void Octree::get_bounding_boxes(BoundingBoxes &boxes) {
	boxes.clear();
    get_bounding_boxes_node(root, boxes);
}

void Octree::get_bounding_boxes_node(OctreeNode *node, BoundingBoxes &boxes) {
	if (node == NULL) return;

	boxes.push_back(node->box);

	for (int i = 0; i < 8; ++i)
		get_bounding_boxes_node(node->child[i], boxes);
}

} // end namespace
