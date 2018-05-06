#ifndef _LIGHTSLICE_H_
#define _LIGHTSLICE_H_

#include "cluster.h"
#include "kdtree_n.h"
#include "light_transport_matrix.h"
#include "kdtree_2.h"
#include "image.h"

#include "events.h"

namespace Renzoku {

class IVirtualPointLightEvaluator;

struct BrdfSample {
    Vec3 w;
    Float pdf;
    Rgb brdf;

    BrdfSample() : pdf(0.0f), w(Vec3(0.0f)) {}
    BrdfSample(const Vec3 &w, Float pdf, const Rgb &brdf) : w(w), pdf(pdf), brdf(brdf) {}
};


struct Slice {
    int start, end;             // start and ending node data
    int representative;
    Rgb color;

    DiscretePdf cluster_pdf;    // how to importance sample a cluster

    ImageFloat luminance_map;
    DiscretePdf2D luminance_map_pdf;

    ImageFloat histogram;       // to record the Markov chain for visualization only

    KdTree2 *kdtree;            // for query nearest neighbors
    NodeData2 *tree_node_data;
    
    vector<BrdfSample> samples;
    int sample_remains;

    //Receiver receiver;          // for debug

    Slice() : kdtree(NULL), tree_node_data(NULL), sample_remains(0) {
    }

    ~Slice() {
        if (tree_node_data) delete[] tree_node_data;
        if (kdtree) delete kdtree;        
    }
};

/**
 * A virtual point light clustering approach based on the approach in 
 * 
 * LightSlice, J. Ou and F. Pellacini, SIGGRAPH Asia 2009.
 *
 * This approach is an extension from Matrix row-column sampling.
 */ 
class LightSlice : public IMouseObserver, public DirectionSampler {
public:
    LightSlice();
    void initialize(Scene *scene, IVirtualPointLightEvaluator *evaluator, BrdfPointLights &vpls);
    
    /**
     * Set number of clusters relative to total VPLs.
     */
    inline void set_num_clusters(int num_clusters);

    inline void set_max_slice_size(int max_slice_size);
    inline void set_num_neighbor_slices(int num_neighbor_slices);
    inline void set_max_clusters_per_slice(int max_slice_clusters);
    
    /**
     * Compute radiance due to a few clusters (out of all clusters found in row-column sampling).
     */
    Rgb gather_clusters(const Receiver &r, IVirtualPointLightEvaluator *vpl);

    /**
     * Importance sample a direction using incident radiance approximated by clusters.
     */
    virtual Rgb sample(Random &rd, const Receiver &r, Vec3 &wi, Float &pdf);
    virtual Float pdf(const Receiver &r, const Vec3 &wi);
    virtual bool is_valid();
    virtual bool is_singular();

private:
    
    /**
     * Build the kd-tree and find slices for all pixels.
     */
    void generate_slices();

private:
    struct IncidentLight {
        Vec3 wi;
        Rgb radiance;
        int cluster_index;
        Float pdf_wi;
        Rgb brdf;

        IncidentLight() : cluster_index(-1), pdf_wi(0.0f) {}
        IncidentLight(const Vec3 &wi, const Rgb &radiance, int cluster_index) : wi(wi), radiance(radiance), cluster_index(cluster_index) {}
    };

    vector<vector<IncidentLight>> global_incident_cache;

private:
    /**
     * Interpolate incident radiance using simple cone bound
     */    
    Float interpolate_radiance_cone(vector<IncidentLight> &cache, int slice, const Vec3 &wi, bool &is_reliable);

    /**
     * Interpolate incident radiance using least square
     */
    Float interpolate_radiance_least_square(const Receiver &r, const Vec3 &wi, int d, const vector<Float> &coeffs, bool &is_reliable);
    void fit_polynomial(const Receiver &r, int degree, vector<IncidentLight> &cache, vector<Float> &coeffs);

private:
    bool use_poly_interpolation;
    int poly_deg;
    vector<Float> coeffs;

private:
    /**
     * Metropolis sampling
     */
    struct SamplingRecord {
        Vec3 w;
        Float f;
        Float sum_f;                                // estimate integration over f using Monte Carlo
        bool start;                                 // mark to generate a start-up sample
        Onb uvn;
        bool is_valid;                              // false when no incoming radiance found
        Float sum_vpl;                              // estimate integration over f using VPLs only
    };
    vector<SamplingRecord> sampling_records;        // cache sampling states at each slice

    int num_sum_samples;                            // number of samples used to estimte sum of f

private:
    Float density_radius;

    Float density_epsilon_ratio;                    // clamp the incoming radiance relative to the sum of incoming radiance from VPL
                                                    // so that the probability is not too small (which causes high variance)

    Float distance_epsilon;                         // ensure slice pdf is not too large.

public:
    void set_density_radius(Float radius);

private:
    struct Mutation {
        enum Type {
            UNIFORM,
            VPL,
            CONE_PERTURB,
            NUM_MUTATIONS
        };
    };
    DiscretePdf mutation_pdf;

    Float perturb_cone_half_angle;
    Float perturb_cone_pdf;
    Float perturb_cos_theta_max;

private:
    Rgb multi_slice_sample(const Receiver &r, Vec3 &wi, Float &pdf_wi);
    Float multi_slice_pdf(const Receiver &r, const Vec3 &wi);

    void metropolis_fetch_cached_sample(int slice, Vec3 &wi, Float &pdf_wi);
    void metropolis_sample(int slice, Vec3 &wi, Float &pdf_wi);
    Float metropolis_get_pdf(int slice, const Vec3 &wi);    

    struct LocalSampler {
        enum Type {
            METROPOLIS, 
            JENSEN
        };
    };
    LocalSampler::Type sampling_type;

    struct LocalDensity {
        enum Type {
            FIXED_RADIUS,
            NEAREST_NEIGHBOR
        };
    };
    LocalDensity::Type radius_type;
    
    Float luminance_map_epsilon;    // for Jensen's method     
    Float metropolis_luminance_epsilon;  // for our method

private:
    int             metropolis_neighbor_slices;
    int             metropolis_neighbor_slices_when_failed;

    vector<int>     metropolis_tmp_slice_indices;             // temporary storage to avoid storage allocation
    vector<Float>   metropolis_tmp_slice_pdf;
    DiscretePdf     metroplis_tmp_discrete_pdf;

private:
    Float estimate_sum(int slice);
    Float estimate_sum_random(Random &rd, int slice, int num_sum_samples);
    Float eval_incoming_radiance(int slice, Vec2 &x);

private:
    Rgb sample_cluster_cone(const Receiver &r, Vec3 &wi, Float &pdf_wi);
    Float pdf_cluster_cone(const Receiver &r, const Vec3 &wi);

    void sample_slice_cluster_cone(int slice, Vec3 &wi, Float &pdf_wi);
    Float pdf_slice_cluster_cone(int slice, const Vec3 &wi);

private:
    Rgb sample_uniform(const Receiver &r, Vec3 &wi, Float &pdf_wi);
    Float pdf_uniform(const Receiver &r, const Vec3 &wi);

private:
    void build_map(int slice);
    void sample_map(int slice, Vec3 &wi, Float &pdf_wi);
    Float pdf_map(int slice, const Vec3 &wi);

    //Rgb jensen_sample(const Receiver &r, Vec3 &wi, Float &pdf_wi);
    //Float jensen_pdf(const Receiver &r, const Vec3 &wi);

    int luminance_map_resolution;

private:
    bool same_local_wi;
    Float view_similarity_cosine_max;

    bool is_view_similar(const Onb &uvn1, const Onb &uvn2);

public:
    Rgb get_slice_color(const Receiver &r);
    Rgb get_slice_sum(const Receiver &r);       // Metropolis normalization factor

private:
    /*
     * Debug functions when mouse is clicked
     */
    void slice_info(int row, int col);
    void build_map_receiver(int row, int col);

private:
    void save_map(int slice);
    void write_to_histogram(int slice, Vec3 &wi);
    void save_histogram(int slice);

    int frame_index;
public:
    void on_frame_end();

public:
    virtual void on_click(MouseButton button, MouseState state, int x, int y);

/**
 * Data for LightSlice clustering
 */
protected:
    Scene *scene;
    Clusters mrcs_clusters;
    BrdfPointLights *all_vpls;

    int num_clusters;

    vector<Slice> slices;
    vector<Clusters> matrix_clusters;

    KdTree6 *kdtree;
    NodeData6 *node_data;                       // spatial and directional coordinates to cluster surfaces
    Float scene_bb_diag;

    int max_slice_size;                         // number of gather points per slice
    int num_neighbor_slices;                    // number of slices to look up to refine global light cluster
    int max_slice_clusters;                     // number of light clusters per slice (including global)

    int num_brdf_samples;                       // importance sampling using lightslice data    
    vector<Float> incident_radiances;
    vector<bool> is_incident_radiance_reliable;
    Float uniform_cone_half_angle;

    vector<BrdfSample> brdf_samples;
    LightTransportMatrix R;

    //MatrixC<Rgb> matI;
    MatrixC<Float> matV;
    IVirtualPointLightEvaluator *evaluator;
};

inline void LightSlice::set_num_clusters(int num_clusters) {
    this->num_clusters = num_clusters;
}

inline void LightSlice::set_max_slice_size(int max_slice_size) {
    this->max_slice_size = max_slice_size;
}

inline void LightSlice::set_num_neighbor_slices(int num_neighbor_slices) {
    this->num_neighbor_slices = num_neighbor_slices;
}

inline void LightSlice::set_max_clusters_per_slice(int max_slice_clusters) {
    if (max_slice_clusters < num_clusters) {
        max_slice_clusters = 2 * num_clusters;
        Log::info() << "Per slice clusters too few. Auto set to 200% of global clusters: " << max_slice_clusters << endn;        
    }
    this->max_slice_clusters = max_slice_clusters;
}

inline void LightSlice::set_density_radius(Float radius) {
    this->density_radius = radius;
}

} // end namespace

#endif