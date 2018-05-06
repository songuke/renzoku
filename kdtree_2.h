#ifndef _KD_TREE_2_H_
#define _KD_TREE_2_H_

#include "common.h"
#include "scene.h"
#include "boundingbox_n.h"
#include "log.h"

#include <queue>
using namespace std;

namespace Renzoku {
    
/**
    * A specialization of KdTreeN class for 2D use. No clustering is supported.
    */

struct KdTreeNode2 {
    int left, right;
    Float plane;			    // location of the splitting plane
    char axis;				    // splitting plane type: 0: X-axis. 1: Y-axis. 2: Z-axis, and so on.

    KdTreeNode2() {
        left = -1;
        right = -1;
        axis = 127;
    }
};

struct NodeData2 {
    Vec2 p;                     // k-dimensional data    
    int index;                  // user data: store the index to some rendering data, e.g., photon, VPLs.
};

/**
* Functor to compare a NodeData with a query point
*/
struct NodeData2Less {
    NodeData2Less(const Vec2 &p) : p(p) {}

    bool operator()(NodeData2 *a, NodeData2 *b) {
        return (a->p - p).squared_length() < (b->p - p).squared_length();
    }

    Vec2 p;     // the query point
};

struct NodeData2Heap
    : priority_queue< NodeData2*, vector<NodeData2 *>, NodeData2Less > {

    NodeData2Heap(const Vec2 &p)
        : priority_queue< NodeData2*, vector<NodeData2 *>, NodeData2Less >(NodeData2Less(p))

    {
    }

};
typedef vector<NodeData2 *> NodeData2Ptrs;

class KdTree2 {
public:
    /**
    * Kd-tree full construction for nearest neighbor search.
    */
    KdTree2(NodeData2 *node_data, int num_nodes, Scene *scene);

    ~KdTree2();

    void find_nearest(const Vec2 &p, int k, NodeData2Heap &nearests) const;
    void find_nearest(const Vec2 &p, Float radius, NodeData2Heap &nearests) const;

private:
    const int dim;

    /**
    * In and out data are ping-pong buffers to store temporary median search result.
    */
    void build_kdtree(int index,
        NodeData2 *in_data, int start, int end, int depth,
        Random &rd, NodeData2 *out_data);

    void find_nearest_rec(int index, const Vec2 &p, int k, NodeData2Heap &heap) const;
    void find_nearest_rec_radius(int index, const Vec2 &p, Float radius_square, NodeData2Heap &heap) const;

private:
    KdTreeNode2 *nodes;
    int num_nodes;
    int max_nodes;
    int *nodes_data_index;              // index for looking up the data array for each kd-tree node

    NodeData2 *nodes_data;              // median-sorted data

private:
    void allocate_nodes(int leaves);    // pre-allocate a set of tree nodes based on the number of leaves
    int get_new_node();                 // fetch a new node from the pool
};

}  // end namespace

#endif
