#include "kdtree.h"

#include "morton.h"
using namespace Nvidia;

#include "scene.h"
#include "log.h"

#include "light_tree.h"

#include <algorithm>
using namespace std;

namespace Renzoku {

static bool less_x(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().x() < p2.org().x();
}

static bool less_y(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().y() < p2.org().y();
}

static bool less_z(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().z() < p2.org().z();
}

static bool equal_x(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().x() == p2.org().x();
}

static bool equal_y(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().y() == p2.org().y();
}

static bool equal_z(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().z() == p2.org().z();
}

static bool greater_x(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().x() > p2.org().x();
}

static bool greater_y(const LightParticle& p1, const LightParticle& p2) {
    return p1.org().y() > p2.org().y();
}

static bool greater_z(const LightParticle& p1, const  LightParticle& p2) {
    return p1.org().z() > p2.org().z();
}

/**
 * O(n) (on average) k-smallest number selection. 
 * 
 * The re-organized array (with correct k-smallest number in place) is stored in out.
 */
static void quick_select_k(Random &rd,
                           bool(*less)(const LightParticle &a, const LightParticle &b),
                           bool(*equal)(const LightParticle &a, const LightParticle &b),
                           bool(*greater)(const LightParticle &a, const LightParticle &b),
                           int k,
                           LightParticles &in, int start, int end, LightParticles &out) {
    if (start > end) return;

    int p = start + (int)(rd() * (end - start + 1));    // must cast to int before +
    
    LightParticle &v = in[p];
    int pos = start;
    for (int i = start; i <= end; ++i) {
        if (less(in[i], v)) {
            out[pos] = in[i];
            pos++;
        }
    }

    int equal_start = pos;    
    for (int i = start; i <= end; ++i) {
        if (equal(in[i], v)) {
            out[pos] = in[i];
            pos++;
        }
    }
    int equal_end = pos - 1;    
    for (int i = start; i <= end; ++i) {
        if (greater(in[i], v)) {
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
        quick_select_k(rd, less, equal, greater, k, out, start, equal_start - 1, in);
    } else {
        for (int i = start; i <= equal_end; ++i)
            in[i] = out[i];
        quick_select_k(rd, less, equal, greater, k, out, equal_end + 1, end, in);
    }
}

/**
 * start, end: inclusive.
 */
void KdTree::build_kdtree(KdTreeNode *node, LightParticles &particles, int start, int end, int depth) {
    if (start > end) return;

    // splitting plane type
    node->type = depth % 3;
    
    // splitting plane location
    int median = (end + start) / 2;
    if (node->type == 0) {
        std::sort(particles.begin() + start, particles.begin() + end + 1, less_x); // [first, last)                
        node->plane = particles[median].org().x();
    } else if (node->type == 1) {
        std::sort(particles.begin() + start, particles.begin() + end + 1, less_y);
        node->plane = particles[median].org().y();
    } else {
        std::sort(particles.begin() + start, particles.begin() + end + 1, less_z);
        node->plane = particles[median].org().z();
    }
    node->particle = particles[median];
    
    if (start <= median - 1) {
        node->left  = get_new_node();       
        build_kdtree(node->left,  particles, start,      median - 1, depth + 1);
    }
    if (median + 1 <= end) {
        node->right = get_new_node();
        build_kdtree(node->right, particles, median + 1, end,        depth + 1);
    }
}


void KdTree::build_kdtree2(KdTreeNode *node, LightParticles &particles, int start, int end, int depth, Random &rd, LightParticles &median_particles) {
    if (start > end) return;

    // splitting plane type
    node->type = depth % 3;
    
    // splitting plane location
    int median = (end + start) / 2;
    if (node->type == 0) {
        quick_select_k(rd, less_x, equal_x, greater_x, median, particles, start, end, median_particles);
        node->plane = median_particles[median].org().x();
    } else if (node->type == 1) {
        quick_select_k(rd, less_y, equal_y, greater_y, median, particles, start, end, median_particles);
        node->plane = median_particles[median].org().y();
    } else {
        quick_select_k(rd, less_z, equal_z, greater_z, median, particles, start, end, median_particles);
        node->plane = median_particles[median].org().z();
    }
    node->particle = median_particles[median];
    
    if (start <= median - 1) {
        node->left  = get_new_node();       
        build_kdtree2(node->left,  median_particles, start,      median - 1, depth + 1, rd, particles);
    }
    if (median + 1 <= end) {
        node->right = get_new_node();
        build_kdtree2(node->right, median_particles, median + 1, end,        depth + 1, rd, particles);
    }
}


static int count_nodes(KdTreeNode *node) {
    if (!node) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

KdTree::KdTree(Photons &photons) : nodes(NULL) {    
    particles.resize(photons.size());
    for (int i = 0; i < photons.size(); ++i) {
        particles[i] = LightParticle(photons[i]);
    }

    allocate_nodes(photons.size());

    if (particles.size() <= 0) {
        root = NULL;        
    } else {
        root = get_new_node();
        build_kdtree(root, particles, 0, photons.size() - 1, 0);
    }
}

KdTree::KdTree(BrdfPointLights &lights) : nodes(NULL) {    
    particles.resize(lights.size());
    for (int i = 0; i < lights.size(); ++i) {
        particles[i] = LightParticle(&lights[i]);
    }

    allocate_nodes(lights.size());

    if (particles.size() <= 0) {
        root = NULL;        
    } else {
        root = get_new_node();
        build_kdtree(root, particles, 0, lights.size() - 1, 0);
    }
}

KdTree::KdTree(BrdfPointLights &lights, Random &rd) : nodes(NULL) {
    particles.resize(lights.size());
    for (int i = 0; i < lights.size(); ++i) {
        particles[i] = LightParticle(&lights[i]);
    }

    allocate_nodes(lights.size());

    if (particles.size() <= 0) {
        root = NULL;        
    } else {
        root = get_new_node();
        median_particles = particles;
        build_kdtree2(root, particles, 0, lights.size() - 1, 0, rd, median_particles);
    }
}

KdTree::KdTree(DensityPoints &points) : nodes(NULL) {    
    particles.resize(points.size());
    for (int i = 0; i < points.size(); ++i) {
        particles[i] = LightParticle(points[i]);
    }

    allocate_nodes(points.size());

    if (particles.size() <= 0) {
        root = NULL;        
    } else {
        root = get_new_node();
        build_kdtree(root, particles, 0, points.size() - 1, 0);
    }
}

KdTree::KdTree(PathNodePtrs &nodes, Scene *scene) : nodes(NULL) {
    particles.resize(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        particles[i] = LightParticle(nodes[i]);
    }
        
    //sort_particles(scene);    // not faster than no sort

    if (particles.size() <= 0) {
        root = NULL;        
    } else {        
        allocate_nodes(nodes.size());
        root = get_new_node();
        median_particles = particles;
        Random &rd = *scene->get_random();
        build_kdtree2(root, particles, 0, nodes.size() - 1, 0, rd, median_particles);
    }
}

KdTree::KdTree(LightNode *light_nodes, int num_nodes, Scene *scene) : nodes(NULL) {
    particles.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        particles[i] = LightParticle(light_nodes + i);
    }

    if (particles.size() <= 0) {
        root = NULL;
    }
    else {
        allocate_nodes(num_nodes);
        root = get_new_node();
        median_particles = particles;
        Random &rd = *scene->get_random();
        build_kdtree2(root, particles, 0, num_nodes - 1, 0, rd, median_particles);
    }
}

void KdTree::allocate_nodes(int leaves) {
    max_nodes = 2 * leaves;
    num_nodes = 0;
    nodes = new KdTreeNode[max_nodes];
}

KdTreeNode *KdTree::get_new_node() {
    if (num_nodes >= max_nodes) {
        Log::info() << "Kd-tree ran out of pre-allocated nodes." << endn;
        return NULL;
    }
    KdTreeNode *node = &nodes[num_nodes++];
    return node;
}

KdTree::~KdTree() {
    if (nodes)
        delete [] nodes;
    root = NULL;
}

static void find_nearest_rec(KdTreeNode *node, const Vec3 &p, int k, LightParticleHeap &heap) {
    if (!node) return;

    if (heap.size() < k) {        
        heap.push(&node->particle);
    } else {
        LightParticle *furthest = heap.top();       
        if ((node->particle.org() - p).squared_length() < (furthest->org() - p).squared_length()) {
            heap.pop();
            heap.push(&node->particle);
        }
    }

    // Find the distance from the photon to the splitting plane
    Float dist; 
    if (node->type == 0) 
        dist = p.x() - node->plane;
    else if (node->type == 1)   
        dist = p.y() - node->plane;
    else 
        dist = p.z() - node->plane;

    // Recur into left and right child node if needed    
    if (dist < 0) { // left
        find_nearest_rec(node->left, p, k, heap);

        LightParticle *furthest = heap.top();
        if (heap.size() < k || (furthest->org() - p).squared_length() > dist * dist) { 
            find_nearest_rec(node->right, p, k, heap);  // check more if not enough, or because the furthest might be replaced by a nearer particle
        }
    } else {        // right
         find_nearest_rec(node->right, p, k, heap);
         
         LightParticle *furthest = heap.top();
         if (heap.size() < k || (furthest->org() - p).squared_length() > dist * dist) {
            find_nearest_rec(node->left, p, k, heap);
         }
    }    
}

void KdTree::find_nearest(const Vec3 &p, int k, LightParticleHeap &nearest_photons) const {    
    find_nearest_rec(root, p, k, nearest_photons);
}

static void find_nearest_rec_radius(KdTreeNode *node, const Vec3 &p, Float radius_square, LightParticleHeap &heap) {
    if (!node) return;

    if ((node->particle.org() - p).squared_length() < radius_square) {        
        heap.push(&node->particle);
    }
    
    // Find the distance from the photon to the splitting plane
    Float dist; 
    if (node->type == 0) 
        dist = p.x() - node->plane;
    else if (node->type == 1)   
        dist = p.y() - node->plane;
    else 
        dist = p.z() - node->plane;

    // Recur into left and right child node if needed    
    if (dist < 0) { // left        
        find_nearest_rec_radius(node->left, p, radius_square, heap);

        if (dist * dist < radius_square)
            find_nearest_rec_radius(node->right, p, radius_square, heap);

    } else {        // right
         find_nearest_rec_radius(node->right, p, radius_square, heap);
         
         if (dist * dist < radius_square)
             find_nearest_rec_radius(node->left, p, radius_square, heap);
    }    
}

void KdTree::find_nearest(const Vec3 &p, Float radius, LightParticleHeap &nearest_particles) const {
    find_nearest_rec_radius(root, p, radius * radius, nearest_particles);
}

void KdTree::sort_particles(Scene *scene) {
    unsigned int *morton = new unsigned int[particles.size()];    
    for (int i = 0; i < particles.size(); ++i) {        
        morton[i] = morton_point(particles[i].org(), &scene->get_bounding_box());
    }
    MortonEntry *entries;
    sort_morton_array(morton, particles.size(), entries);
    LightParticles old_particles;
    old_particles.insert(old_particles.end(), particles.begin(), particles.end());
    for (int i = 0; i < particles.size(); ++i) {
        particles[i] = old_particles[entries[i].index];
    }
    delete [] morton;
    delete [] entries;
}

}  // end namespace
