#ifndef _CLUSTER_BUILDER_H_
#define _CLUSTER_BUILDER_H_

#include "cluster.h"
#include "light_transport_matrix.h"

namespace Renzoku {

/**
 * Build clusters based on the reduced light transport matrix.
 */
class ClusterBuilder {
public:
    ClusterBuilder();
    
    /**     
     * Hasan's approach in SG07.
     */
    void cluster_by_sampling(Scene *scene, 
                            IVirtualPointLightEvaluator *evaluator, 
                            LightTransportMatrix &R, 
                            int num_clusters, Clusters &clusters,                            
                            const BrdfPointLights &all_vpls, 
                            const Pixels &pixels);
    void cluster_by_splitting(Clusters &clusters);

    
    /**
     * Visibility clustering
     *
     * Davidovic's approach in SA10.
     */
    void cluster_visibility(Scene *scene, BrdfPointLights &vpls);

    /**
     * Form factor & brdf clustering
     *
     * Our experimental approach.
     */
    void cluster_form_factor(Scene *scene, BrdfPointLights &vpls);
    void cluster_brdf_light(Scene *scene, BrdfPointLights &vpls);
    void cluster_brdf_surface(Scene *scene, BrdfPointLights &vpls);

};

} // end namespace

#endif