#ifndef _KD_TREE_N_H_
#define _KD_TREE_N_H_

#include "common.h"
#include "scene.h"
#include "boundingbox_n.h"
#include "log.h"

#include <queue>
using namespace std;

namespace Renzoku {

struct KdTreeNodeN {
	int left, right;
	Float plane;			    // location of the splitting plane
	char axis;				    // splitting plane type: 0: X-axis. 1: Y-axis. 2: Z-axis, and so on.
    
    bool clustered;             // is this the node where clustering is done?

	KdTreeNodeN() {
		left = -1;
		right = -1;
        plane = FLT_MAX;
		axis = 127;	
        clustered = false;	
	}
};

template <int dim>
struct NodeData {
    // input
    VecN<dim> p;                // k-dimensional data
    int index;                  // store the index to some rendering data, e.g., photon, VPLs.
                                
    // output
    int cluster;                // when the kd-tree is used for clustering, this will be used to store cluster index.

    NodeData() {
        cluster = -1;
    }
};

template <int dim>
struct ClusterData {    
    int start, end;             // index of starting and ending node data in the cluster
    BoundingBoxN<dim> bb;       // data bounding box

    ClusterData(int start, int end, const BoundingBoxN<dim> &bb) : start(start), end(end), bb(bb) {
    }
};

/**
 * Functor to compare a NodeData with a query point
 */
template <int dim>
struct NodeDataLess {
    NodeDataLess(const VecN<dim> &p) : p(p) {}

    bool operator()(NodeData<dim> *a, NodeData<dim> *b) {
        return (a->p - p).squared_length() < (b->p - p).squared_length();
    }

    VecN<dim> p;     // the query point
};

template <int dim>
struct NodeDataHeap 
        : priority_queue< NodeData<dim>*, vector<NodeData<dim> *>, NodeDataLess<dim> > {

    NodeDataHeap(const VecN<dim> &p) 
        : priority_queue< NodeData<dim>*, vector<NodeData<dim> *>, NodeDataLess<dim> >( NodeDataLess<dim>(p) ) 
    
    {
    }

};

template <int dim>
struct ClusterIndexLess {
    ClusterIndexLess(const vector<ClusterData<dim>> &clusters, const VecN<dim> &p) : clusters(clusters), p(p) {        
    }

    bool operator()(int a, int b) const {
        return clusters[a].bb.min_squared_distance(p) < clusters[b].bb.min_squared_distance(p);
    }

    const vector<ClusterData<dim>> &clusters;
    VecN<dim> p;     // the query point
};

template <int dim>
struct ClusterIndexHeap 
        : priority_queue< int, vector<int>, ClusterIndexLess<dim> > {

    ClusterIndexHeap(const vector<ClusterData<dim>> &clusters, const VecN<dim> &p) 
        : priority_queue< int, vector<int>, ClusterIndexLess<dim> >( ClusterIndexLess<dim>(clusters, p) ) 
    
    {
    }
};

template <int dim>
class KdTreeN {
public:
	/**
     * Kd-tree full construction for nearest neighbor search.
     */
    KdTreeN(NodeData<dim> *node_data, int num_nodes, Scene *scene);
    
    /**
     * Kd-tree partial construction for clustering. 
     * We build the tree top-down, but stop split a node when the number of nodes is less than a threshold.
     *
     * All data that belong to a node becomes a cluster.
     */
    KdTreeN(int max_clusters, 
            NodeData<dim> *node_data, int num_nodes, Scene *scene,
            vector<pair<int, int>> &clusters);

	~KdTreeN();	
        
    void find_nearest(const VecN<dim> &p, int k, NodeDataHeap<dim> &nearests) const;
    void find_nearest(const VecN<dim> &p, Float radius, NodeDataHeap<dim> &nearests) const;
    
    void find_cluster(const VecN<dim> &p, int &cluster_index) const;

    /**
     * Furthest cluster is at position 0.
     */
    void find_nearest_cluster(const VecN<dim> &p, int k, vector<int> &cluster_index) const;

    void get_cluster_bounding_boxes(BoundingBoxes &bbs) const;

private:
    /**
     * In and out data are ping-pong buffers to store temporary median search result.
     */
    void build_kdtree(int index, 
                      NodeData<dim> *in_data, int start, int end, int depth, 
                      Random &rd, NodeData<dim> *out_data);

    void cluster_kdtree(int max_cluster_size, vector<pair<int, int>> &clusters,
                        int index, NodeData<dim> *in_data, int start, int end, int depth, 
                        Random &rd, NodeData<dim> *out_data); 

    void find_nearest_rec(int index, const VecN<dim> &p, int k, NodeDataHeap<dim> &heap) const;
    void find_nearest_rec_radius(int index, const VecN<dim> &p, Float radius_square, NodeDataHeap<dim> &heap) const;

    void find_cluster_rec(int index, const VecN<dim> &p, int &cluster_index) const;
    void find_nearest_cluster_rec(int index, const VecN<dim> &p, int k, ClusterIndexHeap<dim> &cluster_index) const;

private:
	KdTreeNodeN *nodes;
    int num_nodes;
    int max_nodes;
    int *nodes_data_index;              // index for looking up the data array for each kd-tree node
        
    NodeData<dim> *nodes_data;          // median-sorted data

    vector<ClusterData<dim>> cluster_bbs;    // cluster centroid for nearest cluster search

    static const Float PLANE_BIAS;

private:
    void allocate_nodes(int leaves);    // pre-allocate a set of tree nodes based on the number of leaves
    int get_new_node();                 // fetch a new node from the pool
};

template <int dim>
const Float KdTreeN<dim>::PLANE_BIAS = -1e-3f;      // shift the plane a little to avoid numerical inaccuracy 
                                                    // when testing a point at the left or right of the splitting plane)
                                                    // this inaccuracy can appear as noise in surface clustering, e.g., in LightSlice,
                                                    // (display surface cluster as color to see)
                                                    // 
                                                    // we only use this when we need to determine the cluster for point
                                            


template <int dim>
static int compare_data(int axis, const NodeData<dim> &a, const NodeData<dim> &b) {
    return (a.p[axis] < b.p[axis]) ? -1 : 
           (a.p[axis] == b.p[axis]) ? 0 : 1;
}

/**
 * O(n) (on average) k-smallest number selection. 
 * 
 * The re-organized array (with correct k-smallest number in place) is stored both in in/out at the end.
 */
template <int dim>
static void quick_select_k(Random &rd,
                           int axis,
                           int (*compare)(int axis, const NodeData<dim> &a, const NodeData<dim> &b),
                           int k,
                           NodeData<dim> *in, int start, int end, NodeData<dim> *out) {
    if (start > end) return;

    int p = start + (int)(rd() * (end - start + 1));    // must cast to int before +
    
    NodeData<dim> &v = in[p];
    int pos = start;
    // FIXME: we can do an in-place swap and discard the out-array
    for (int i = start; i <= end; ++i) {
        if (compare(axis, in[i], v) < 0) {
            out[pos] = in[i];
            pos++;
        }
    }

    int equal_start = pos;    
    for (int i = start; i <= end; ++i) {
        if (compare(axis, in[i], v) == 0) {
            out[pos] = in[i];
            pos++;
        }
    }
    int equal_end = pos - 1;    
    for (int i = start; i <= end; ++i) {
        if (compare(axis, in[i], v) > 0) {
            out[pos] = in[i];
            pos++;
        }
    }

    // ignoring the entire equal part is important to make the search fast
    if (k >= equal_start && k <= equal_end) {
        // sync back the unused part
        for (int i = start; i <= end; ++i)
            in[i] = out[i];     
        return;
    } else if (k < equal_start) {
        for (int i = equal_start; i <= end; ++i)
            in[i] = out[i];
        quick_select_k(rd, axis, compare, k, out, start, equal_start - 1, in);
    } else {
        for (int i = start; i <= equal_end; ++i)
            in[i] = out[i];
        quick_select_k(rd, axis, compare, k, out, equal_end + 1, end, in);
    }
}

template <int dim>
void KdTreeN<dim>::build_kdtree(int index, NodeData<dim> *data, int start, int end, int depth, Random &rd, NodeData<dim> *median_data) {
        
    //if (start > end) return;
    
    // Apr 05: 
    if (start >= end) return;
    
    KdTreeNodeN *node = &nodes[index];

    // take node's bounding box and determine the longest axis
    BoundingBoxN<dim> bb;
    for (int i = start; i <= end; ++i)
        bb.merge(data[i].p);
    node->axis = bb.get_longest_axis();

    // splitting plane location
    int median = (end + start) / 2;

    quick_select_k<dim>(rd, node->axis, compare_data<dim>, median, data, start, end, median_data);
    node->plane = median_data[median].p[node->axis];
    nodes_data_index[index] = median;

    // Apr 05: must include median
    if (start <= median) {

    //if (start <= median - 1) {
    
        node->left  = get_new_node();

        //build_kdtree(node->left,  median_data, start,      median - 1, depth + 1, rd, data);

        // Apr 05: must include median        
        build_kdtree(node->left,  median_data, start,      median, depth + 1, rd, data);
    }
    if (median + 1 <= end) {
        node->right = get_new_node();        
        build_kdtree(node->right, median_data, median + 1, end,        depth + 1, rd, data);
    }
}

template <int dim>
void KdTreeN<dim>::cluster_kdtree(
                      int max_cluster_size, vector<pair<int, int>> &clusters, 
                      int index, 
                      NodeData<dim> *data, int start, int end, 
                      int depth, Random &rd, 
                      NodeData<dim> *median_data) {

    // Apr 05 fix
    if (start >= end) return;
    
    //if (start > end) return;    

    KdTreeNodeN *node = &nodes[index];

    // Capture cluster at some level, and then continue to build the full-tree.
    // It is possible to stop the kd-tree build at this level,
    // but the tree is useless for accurate nearest neighbor search.
    if (end - start + 1 <= max_cluster_size && 
        data[start].cluster < 0 && 
        node->clustered == false) {      
                          
        // store cluster index into the NodeData
        for (int i = start; i <= end; ++i) {

            // just sync the output data because the node data is updated        
            data[i].cluster = clusters.size();            
            median_data[i].cluster = clusters.size();
        }

        // cluster data to return to user
        clusters.push_back(pair<int, int>(start, end));
        
        // internally cache the bounding box of the cluster for nearest cluster search
        BoundingBoxN<dim> bb;
        for (int i = start; i <= end; ++i) {            
            bb.merge(data[i].p);
        }
        ClusterData<dim> cluster(start, end, bb);        
        this->cluster_bbs.push_back(cluster);

        node->clustered = true;        
    }
    

    // take node's bounding box and determine the longest axis
    BoundingBoxN<dim> bb;
    for (int i = start; i <= end; ++i)
        bb.merge(data[i].p);
    node->axis = bb.get_longest_axis();

    // splitting plane location
    int median = (end + start) / 2;

    quick_select_k<dim>(rd, node->axis, compare_data<dim>, median, data, start, end, median_data);
    node->plane = median_data[median].p[node->axis];
    nodes_data_index[index] = median;
    
    // Apr 05: must include median
    if (start <= median) {
    
    //if (start <= median - 1) {    
        node->left  = get_new_node();
        nodes[node->left].clustered = node->clustered;    // no more clustering if this node is already clustered

        //cluster_kdtree(max_cluster_size, clusters, node->left,  median_data, start,      median - 1, depth + 1, rd, data);
        
        // Apr 05:
        cluster_kdtree(max_cluster_size, clusters, node->left,  median_data, start,      median, depth + 1, rd, data);
    }
    if (median + 1 <= end) {
        node->right = get_new_node();        
        nodes[node->right].clustered = node->clustered;

        cluster_kdtree(max_cluster_size, clusters, node->right, median_data, median + 1, end,        depth + 1, rd, data);
    }
}

/**
 * Full kd-tree construction for nearest neighbor search.
 *
 * Node data can be modified.
 */
template <int dim>
KdTreeN<dim>::KdTreeN(NodeData<dim> *data, int num_nodes, Scene *scene) : nodes(NULL), nodes_data_index(NULL), nodes_data(NULL) {    
    if (num_nodes > 0) {
        NodeData<dim> *median_data = new NodeData<dim>[num_nodes];
        memcpy(median_data, data, sizeof(NodeData<dim>) * num_nodes);

        allocate_nodes(num_nodes);
        int root = get_new_node();

        Random &rd = *scene->get_random();
        build_kdtree(root, data, 0, num_nodes - 1, 0, rd, median_data);

        // TODO: test:
        // after build, the data and median_data is always in sync.
        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (median_data[i].p[j] != data[i].p[j]) {
                    Log::error() << "Node " << i << ": " << "Data not in sync." << endn;
                }
            }
        }
        
        delete [] median_data;

        this->nodes_data = data;
    }
}

/**
 * Full kd-tree construction with clustering information captured.
 */
template <int dim>
KdTreeN<dim>::KdTreeN(int max_cluster_size, 
                      NodeData<dim> *data, int num_nodes, Scene *scene,
                      vector<pair<int, int>> &clusters) : nodes(NULL), nodes_data_index(NULL), nodes_data(NULL) {    

    if (num_nodes > 0) {

        NodeData<dim> *median_data = new NodeData<dim>[num_nodes];
        memcpy(median_data, data, sizeof(NodeData<dim>) * num_nodes);

        allocate_nodes(num_nodes);
        int root = get_new_node();

        Random &rd = *scene->get_random();
        cluster_kdtree(max_cluster_size, clusters, 
                root, data, 0, num_nodes - 1, 
                0, rd, 
                median_data);

        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (median_data[i].p[j] != data[i].p[j]) {
                    Log::error() << "Node " << i << ": " << "Data not in sync." << endn;
                }
            }
        }
        
        delete [] median_data;

        this->nodes_data = data;
    }
}

template <int dim>
void KdTreeN<dim>::allocate_nodes(int leaves) {
    max_nodes = 2 * leaves;
    num_nodes = 0;
    nodes = new KdTreeNodeN[max_nodes];
    nodes_data_index = new int[max_nodes];
}

template <int dim>
int KdTreeN<dim>::get_new_node() {
    if (num_nodes >= max_nodes) {
        Log::info() << "Kd-tree ran out of pre-allocated nodes." << endn;
        return -1;
    }
    int index = num_nodes;
    num_nodes++;
    return index;
}

template <int dim>
KdTreeN<dim>::~KdTreeN() {
    if (nodes) {
        delete [] nodes;
        delete [] nodes_data_index;
    }
}

template <int dim>
void KdTreeN<dim>::find_nearest_rec(int index, const VecN<dim> &p, int k, NodeDataHeap<dim> &heap) const {
    if (index < 0) return;

    if (heap.size() < k) {        
        heap.push(&nodes_data[nodes_data_index[index]]);
    } else {
        NodeData<dim> *furthest = heap.top();       
        if ((nodes_data[nodes_data_index[index]].p - p).squared_length() < (furthest->p - p).squared_length()) {
            heap.pop();
            heap.push(&nodes_data[nodes_data_index[index]]);
        }
    }
        
    // Find the distance from the photon to the splitting plane
    KdTreeNodeN *node = &nodes[index];

    // no splitting then return
    if (node->axis >= dim) return;

    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane;    
    if (dist <= 0.0f) { // left
        find_nearest_rec(node->left, p, k, heap);

        NodeData<dim> *furthest = heap.top();
        if (heap.size() < k || (furthest->p - p).squared_length() > dist * dist) { 
            find_nearest_rec(node->right, p, k, heap);  // check more if not enough, or because the furthest might be replaced by a nearer particle
        }
    } else {        // right
         find_nearest_rec(node->right, p, k, heap);
         
         NodeData<dim> *furthest = heap.top();
         if (heap.size() < k || (furthest->p - p).squared_length() > dist * dist) {
            find_nearest_rec(node->left, p, k, heap);
         }
    }    
}

template <int dim>
void KdTreeN<dim>::find_nearest(const VecN<dim> &p, int k, NodeDataHeap<dim> &nearests) const {    
    if (num_nodes > 0)
        find_nearest_rec(0, p, k, nearests);
}

template <int dim>
void KdTreeN<dim>::find_nearest_rec_radius(int index, const VecN<dim> &p, Float radius_square, NodeDataHeap<dim> &heap) const {
    if (index < 0) return;

    if ((nodes_data[nodes_data_index[index]].p - p).squared_length() < radius_square) {    
        // FIXME: in fact we don't need a heap here if the output does not need to be in order from furthest to nearest
        heap.push(&nodes_data[nodes_data_index[index]]);
    }
        
    KdTreeNodeN *node = &nodes[index];

    if (node->axis >= dim) return;
    
    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane; 
    if (dist <= 0.0f) { // left        
        find_nearest_rec_radius(node->left, p, radius_square, heap);

        if (dist * dist < radius_square)
            find_nearest_rec_radius(node->right, p, radius_square, heap);

    } else {        // right
         find_nearest_rec_radius(node->right, p, radius_square, heap);
         
         if (dist * dist < radius_square)
             find_nearest_rec_radius(node->left, p, radius_square, heap);
    }    
}

template <int dim>
void KdTreeN<dim>::find_nearest(const VecN<dim> &p, Float radius, NodeDataHeap<dim> &nearests) const {
    if (num_nodes > 0)
        find_nearest_rec_radius(0, p, radius * radius, nearests);
}


template <int dim>
void KdTreeN<dim>::find_cluster_rec(int index, const VecN<dim> &p, int &cluster_index) const {
    
    // Find the distance from the photon to the splitting plane
    KdTreeNodeN *node = &nodes[index];

    // leaf node (cluster node)
    if (node->clustered) {
        cluster_index = nodes_data[nodes_data_index[index]].cluster;
        return;
    }

    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane + PLANE_BIAS;
    
    if (dist <= 0.0f) {
        find_cluster_rec(node->left, p, cluster_index);
    }
    else {

        find_cluster_rec(node->right, p, cluster_index);

    }
}


template <int dim>
void KdTreeN<dim>::find_cluster(const VecN<dim> &p, int &cluster_index) const {
    if (num_nodes > 0)
        find_cluster_rec(0, p, cluster_index);
    else
        cluster_index = -1;
}


template <int dim>
void KdTreeN<dim>::find_nearest_cluster_rec(int index, const VecN<dim> &p, int k, ClusterIndexHeap<dim> &heap) const {
    
    KdTreeNodeN *node = &nodes[index];

    if (node->clustered) {
        int cluster_index = nodes_data[nodes_data_index[index]].cluster;
        if (heap.size() < k) {        
            heap.push(cluster_index);
        } else {
            int furthest_index = heap.top();
            if (cluster_bbs[cluster_index].bb.min_squared_distance(p) < cluster_bbs[furthest_index].bb.min_squared_distance(p)) {
                heap.pop();
                heap.push(cluster_index);
            }
        }

        return;
    }

    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane + PLANE_BIAS;
    if (dist <= 0.0f) {
        find_nearest_cluster_rec(node->left, p, k, heap);
        
        int furthest_index = heap.top();
        if (heap.size() < k || cluster_bbs[furthest_index].bb.min_squared_distance(p) > dist * dist) {
            find_nearest_cluster_rec(node->right, p, k, heap);
        }

    } else {
        find_nearest_cluster_rec(node->right, p, k, heap);
         
        int furthest_index = heap.top();
        if (heap.size() < k || cluster_bbs[furthest_index].bb.min_squared_distance(p) > dist * dist) {
            find_nearest_cluster_rec(node->left, p, k, heap);
        }
    }    
}


template <int dim>
void KdTreeN<dim>::find_nearest_cluster(const VecN<dim> &p, int k, vector<int> &cluster_index) const {

    /*
    // test brute force cluster search
    cluster_index.clear();
    Float min_dist = FLT_MAX;
    int index = -1;
    for (int i = 0; i < cluster_bbs.size(); ++i) {
        Float dist = ((cluster_bbs[i].bb.vmin + cluster_bbs[i].bb.vmax)*0.5f - p).squared_length();
        if (dist < min_dist) {
            min_dist = dist;
            index = i;
        }
    }
    if (index >= 0) {
        cluster_index.push_back(index);
    }
    return;
    */

    cluster_index.clear();
    if (num_nodes > 0) {

        ClusterIndexHeap<dim> heap(cluster_bbs, p);
        find_nearest_cluster_rec(0, p, k, heap);
                
        while (heap.size() > 0) {
            cluster_index.push_back(heap.top());
            heap.pop();
        }
    }
}

template <int dim>
void KdTreeN<dim>::get_cluster_bounding_boxes(BoundingBoxes &bbs) const {
    bbs.clear();
    for (int i = 0; i < cluster_bbs.size(); ++i) {
        BoundingBox bb;
        // only take the first three dimensions
        bb.v_min = Vec3(cluster_bbs[i].bb.vmin[0], cluster_bbs[i].bb.vmin[1], cluster_bbs[i].bb.vmin[2]);
        bb.v_max = Vec3(cluster_bbs[i].bb.vmax[0], cluster_bbs[i].bb.vmax[1], cluster_bbs[i].bb.vmax[2]);
        bbs.push_back(bb);
    }
}

typedef VecN<6> Vec6;
typedef NodeData<6> NodeData6;
typedef KdTreeN<6> KdTree6;
typedef NodeDataHeap<6> NodeData6Heap;

typedef VecN<9> Vec9;
typedef NodeData<9> NodeData9;
typedef KdTreeN<9> KdTree9;
typedef NodeDataHeap<9> NodeData9Heap;

}  // end namespace

#endif
