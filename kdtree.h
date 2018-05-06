#ifndef _KD_TREE_H_
#define _KD_TREE_H_

#include "common.h"
#include "light_particle.h"

#include <queue>
using namespace std;

namespace Renzoku {
    
struct KdTreeNode {
	KdTreeNode *left, *right;
	Float plane;			    // location of the splitting plane
	char type;				    // splitting plane type: 0: X-axis. 1: Y-axis. 2: Z-axis. 

    LightParticle particle;		// the particle array can be in any order during sort, so need to store data as a copy.		
    
	KdTreeNode() {
		left = NULL;
		right = NULL;
		type = 0;		
	}
};

/**
 * Special kd-tree with k = 3. 
 *
 * The plane to split is in turn, x, y, z. 
 *
 * For a general kd-tree that handles n-dimensional data, use KdTreeN class.
 */
class KdTree {
public:
	KdTree(Photons &photons);           
	KdTree(BrdfPointLights &vpls);
    KdTree(DensityPoints &points);    

    /**
     * Random generator is used for finding the median using a randomized algorithm.
     */
    KdTree(BrdfPointLights &vpls, Random &rd);
    KdTree(PathNodePtrs &nodes, Scene *scene);
    KdTree(LightNode *nodes, int num_nodes, Scene *scene);

	~KdTree();	
    
    void find_nearest(const Vec3 &p, int k, LightParticleHeap &nearest_particles) const;
    void find_nearest(const Vec3 &p, Float radius, LightParticleHeap &nearest_particles) const;

private:
    void build_kdtree (KdTreeNode *node, LightParticles &particles, int start, int end, int depth);
    void build_kdtree2(KdTreeNode *node, LightParticles &particles, int start, int end, int depth, 
                       Random &rd, LightParticles &median_particles);
    
    /**
     * Sort points by Morton code
     */
    void sort_particles(Scene *scene);

private:
	KdTreeNode *root;
    LightParticles particles;
    LightParticles median_particles;	// to store linear-time median search result
    
private:
    void allocate_nodes(int leaves);    // pre-allocate a set of tree nodes based on the number of leaves
    KdTreeNode *get_new_node();         // fetch a new node from the pool

    KdTreeNode *nodes;
    int num_nodes;
    int max_nodes;
};

}  // end namespace

#endif
