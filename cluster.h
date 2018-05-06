#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include "pdf.h"
#include "cone.h"
#include "rgb.h"

#include <vector>
using namespace std;

namespace Renzoku {

/**
 * Cluster represents how rows/columns of a light transport matrix can be grouped together.
 * 
 * The index it stores refers to the global row/column index in the matrix.
 */
struct Cluster {
    int center;                 // cluster center
    int representative;         
    Float weight;               // accounts for the probability of choosing the representative
    vector<int> indices;        
    DiscretePdf pdf;            // the probability distribution of all columns in the cluster

    Cone bc;                          // bounding cone of all VPLs to the slice representative
    Rgb total_incident_radiance;      // incoming radiance to a slice, with visibility
    Rgb total_outgoing_radiance;

    Vec3 slice_wi;
    bool is_upper_hemisphere;           // true if the cluster representative is on the upper part of the hemisphere at the slice

    Cluster() : bc(Vec3(0.0f), 0.0f), is_upper_hemisphere(false) {
        center = -1;
        representative = -1;
        weight = 0;    
    }
};

typedef vector<Cluster> Clusters;

} // end namespace

#endif