#include "kdtree_2.h"

namespace Renzoku {

static int compare_data(int axis, const NodeData2 &a, const NodeData2 &b) {
    return (a.p[axis] < b.p[axis]) ? -1 :
        (a.p[axis] == b.p[axis]) ? 0 : 1;
}

/**
* O(n) (on average) k-smallest number selection.
*
* The re-organized array (with correct k-smallest number in place) is stored both in in/out at the end.
*/
static void quick_select_k(Random &rd,
    int axis,
    int(*compare)(int axis, const NodeData2 &a, const NodeData2 &b),
    int k,
    NodeData2 *in, int start, int end, NodeData2 *out) {

    if (start > end) return;

    int p = start + (int)(rd() * (end - start + 1));    // must cast to int before +

    NodeData2 &v = in[p];
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
    }
    else if (k < equal_start) {
        for (int i = equal_start; i <= end; ++i)
            in[i] = out[i];
        quick_select_k(rd, axis, compare, k, out, start, equal_start - 1, in);
    }
    else {
        for (int i = start; i <= equal_end; ++i)
            in[i] = out[i];
        quick_select_k(rd, axis, compare, k, out, equal_end + 1, end, in);
    }
}

void KdTree2::build_kdtree(int index, NodeData2 *data, int start, int end, int depth, Random &rd, NodeData2 *median_data) {
    if (start > end) return;

    KdTreeNode2 *node = &nodes[index];

    // take node's bounding box and determine the longest axis
    BoundingBoxN<2> bb;
    VecN<2> pos;
    for (int i = start; i <= end; ++i) {
        pos[0] = data[i].p.x();
        pos[1] = data[i].p.y();
        bb.merge(pos);
    }
    node->axis = bb.get_longest_axis();

    // splitting plane location
    int median = (end + start) / 2;

    quick_select_k(rd, node->axis, compare_data, median, data, start, end, median_data);
    node->plane = median_data[median].p[node->axis];
    nodes_data_index[index] = median;

    if (start <= median - 1) {
        node->left = get_new_node();
        build_kdtree(node->left, median_data, start, median - 1, depth + 1, rd, data);
    }
    if (median + 1 <= end) {
        node->right = get_new_node();
        build_kdtree(node->right, median_data, median + 1, end, depth + 1, rd, data);
    }
}


/**
* Full kd-tree construction for nearest neighbor search.
*
* Node data can be modified.
*/
KdTree2::KdTree2(NodeData2 *data, int num_nodes, Scene *scene) : dim(2), nodes(NULL), nodes_data_index(NULL), nodes_data(NULL) {
    if (num_nodes > 0) {
        NodeData2 *median_data = new NodeData2[num_nodes];
        memcpy(median_data, data, sizeof(NodeData2) * num_nodes);

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

        delete[] median_data;

        this->nodes_data = data;
    }
}

void KdTree2::allocate_nodes(int leaves) {
    max_nodes = 2 * leaves;
    num_nodes = 0;
    nodes = new KdTreeNode2[max_nodes];
    nodes_data_index = new int[max_nodes];
}

int KdTree2::get_new_node() {
    if (num_nodes >= max_nodes) {
        Log::info() << "Kd-tree ran out of pre-allocated nodes." << endn;
        return -1;
    }
    int index = num_nodes;
    num_nodes++;
    return index;
}

KdTree2::~KdTree2() {
    if (nodes) {
        delete[] nodes;
        delete[] nodes_data_index;
    }
}

void KdTree2::find_nearest_rec(int index, const Vec2 &p, int k, NodeData2Heap &heap) const {
    if (index < 0) return;

    if (heap.size() < k) {        
        heap.push(&nodes_data[nodes_data_index[index]]);
    }
    else {
        NodeData2 *furthest = heap.top();
        if ((nodes_data[nodes_data_index[index]].p - p).squared_length() < (furthest->p - p).squared_length()) {
            heap.pop();
            heap.push(&nodes_data[nodes_data_index[index]]);
        }
    }

    // Find the distance from the photon to the splitting plane
    KdTreeNode2 *node = &nodes[index];

    // no splitting then return
    if (node->axis >= dim) return;

    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane;
    if (dist < 0) { // left
        find_nearest_rec(node->left, p, k, heap);

        NodeData2 *furthest = heap.top();
        if (heap.size() < k || (furthest->p - p).squared_length() > dist * dist) {
            find_nearest_rec(node->right, p, k, heap);  // check more if not enough, or because the furthest might be replaced by a nearer particle
        }
    }
    else {        // right
        find_nearest_rec(node->right, p, k, heap);

        NodeData2 *furthest = heap.top();
        if (heap.size() < k || (furthest->p - p).squared_length() > dist * dist) {
            find_nearest_rec(node->left, p, k, heap);
        }
    }
}

void KdTree2::find_nearest(const Vec2 &p, int k, NodeData2Heap &nearests) const {
    if (num_nodes > 0)
        find_nearest_rec(0, p, k, nearests);
}

void KdTree2::find_nearest_rec_radius(int index, const Vec2 &p, Float radius_square, NodeData2Heap &heap) const {
    if (index < 0) return;

    if ((nodes_data[nodes_data_index[index]].p - p).squared_length() < radius_square) {        
        // FIXME: in fact we don't need a heap here if the output does not need to be in order from furthest to nearest
        heap.push(&nodes_data[nodes_data_index[index]]);
    }

    KdTreeNode2 *node = &nodes[index];

    if (node->axis >= dim) return;

    // Recur into left and right child node if needed    
    Float dist = p[node->axis] - node->plane;
    if (dist < 0) { // left        
        find_nearest_rec_radius(node->left, p, radius_square, heap);

        if (dist * dist < radius_square)
            find_nearest_rec_radius(node->right, p, radius_square, heap);

    }
    else {        // right
        find_nearest_rec_radius(node->right, p, radius_square, heap);

        if (dist * dist < radius_square)
            find_nearest_rec_radius(node->left, p, radius_square, heap);
    }
}

void KdTree2::find_nearest(const Vec2 &p, Float radius, NodeData2Heap &nearests) const {
    if (num_nodes > 0)
        find_nearest_rec_radius(0, p, radius * radius, nearests);
}

}  // end namespace