#ifndef _MATRIX_ROW_COLUMN_SAMPLING_H_
#define _MATRIX_ROW_COLUMN_SAMPLING_H_

#include "cluster.h"

namespace Renzoku {

class IVirtualPointLightEvaluator;

/**
 * A virtual point light clustering approach based on matrix row-column samping
 * approach in 
 * 
 * Matrix Row-Column Sampling for the Many-light Problem, by M. Hasan, F. Pellacini, and K. Bala, SIGGRAPH 2007. 
 */ 
class MatrixRowColumnSampling {
public:
    MatrixRowColumnSampling();
    void initialize(Scene *scene, IVirtualPointLightEvaluator *evaluator, BrdfPointLights &vpls);
    
    inline void set_block_size(int block_size);

    /**
     * Set number of clusters relative to total VPLs.
     */
    inline void set_cluster_ratio(Float ratio);
    inline void set_max_clusters(int max_clusters);

    /**
     * Compute radiance due to a few clusters (out of all clusters found in row-column sampling).
     */
    Rgb gather_clusters(const Receiver &r, IVirtualPointLightEvaluator *vpl);

protected:
    Scene *scene;
    Clusters all_clusters;
    BrdfPointLights *all_vpls;

    int block_size;
    int num_clusters;
    Float cluster_ratio;
    int max_clusters;
};

inline void MatrixRowColumnSampling::set_cluster_ratio(Float ratio) {
    cluster_ratio = ratio;
}

inline void MatrixRowColumnSampling::set_block_size(int block_size) {
    this->block_size = block_size;
}

inline void MatrixRowColumnSampling::set_max_clusters(int max_clusters) {
    this->max_clusters = max_clusters;
}

} // end namespace

#endif