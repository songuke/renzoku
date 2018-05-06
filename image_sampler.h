#ifndef _IMAGE_SAMPLER_H_
#define _IMAGE_SAMPLER_H_

#include "vec2.h"

namespace Renzoku {

typedef vector<Vec2> Pixels;

/**
 * Generate pixel samples that form the reduced light transport matrix
 * which serves the exploration step for column clustering.
 */
class ImageSampler {
public:
    /**
     * Subdivide the image into blocks and for each block pick a pixel.
     *
     * Proposed in Hasan SG07.
     */
    void sample_stratified(Scene *scene, int block_size, Pixels &pixels);

    /**
     * Cluster pixels by building the kd-tree on 6D samples (position, normal) of surface points.
     * 
     * Proposed in Ou SA11. 
     */
    void sample_kdtree(Scene *scene, Pixels &pixels);
};

} // end namespace

#endif
