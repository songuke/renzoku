#include "cluster_builder.h"
#include "stats.h"
#include "random.h"
#include "scene.h"
#include "aggregate.h"
#include "camera.h"
#include "log.h"

namespace Renzoku {

ClusterBuilder::ClusterBuilder() {
}

static Float magnitude(const vector<Float> &v) {
    Float sum = 0.0f;
    for (int i = 0; i < v.size(); ++i)
        sum += v[i] * v[i];
    return sqrtf(sum);
}

/**
 * This only evaluates radiance. 
 *
 * NOTE: It does not support incident radiance clustering at this moment.
 */
/*
static void eval_column(Scene *scene,
                        IVirtualPointLightEvaluator *evaluator, 
                        const BrdfPointLight &vpl, const Pixels &pixels, vector<Float> &c) {
    Camera *camera = scene->get_camera();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Lights &lights = *scene->get_lights();
        
    for (int ridx = 0; ridx < pixels.size(); ++ridx) {
        Vec2 pixel = pixels[ridx];
        
        Ray r = camera->shoot_ray(pixel.y(), pixel.x());
        HitRecord rec;
        if (! agg->hit(r, tmin, max_tmax, tick, rec)) {
            c[3 * ridx    ] = 0.0f;
            c[3 * ridx + 1] = 0.0f;
            c[3 * ridx + 2] = 0.0f;
            continue;
        }
                
        if (rec.light) {
            c[3 * ridx    ] = 0.0f;
            c[3 * ridx + 1] = 0.0f;
            c[3 * ridx + 2] = 0.0f;
            continue;
        }

        Receiver recv(r, rec);                            
        Rgb radiance = evaluator->radiance(scene, recv, vpl);        
        c[3 * ridx    ] = radiance.red();
        c[3 * ridx + 1] = radiance.green();
        c[3 * ridx + 2] = radiance.blue();
    }    
}*/

static Float get_distance(const vector<Float> &mi, const vector<Float> &mj) {    
    Float ci_dot_cj = 0.f;
    for (int r = 0; r < mi.size(); ++r) {        
        ci_dot_cj += mi[r] * mj[r];
    }
    
    return fabs(magnitude(mi) * magnitude(mj) - ci_dot_cj); // do fabs to ensure positive distance (otherwise negative values can appear due to numerical errors).
}

static Float dot(const vector<Float> &mi, const vector<Float> &mj) {    
    Float ci_dot_cj = 0.f;
    for (int r = 0; r < mi.size(); ++r) {        
        ci_dot_cj += mi[r] * mj[r];
    }
    return ci_dot_cj;
}

void ClusterBuilder::cluster_by_sampling(Scene *scene, 
                                         IVirtualPointLightEvaluator *evaluator, 
                                         LightTransportMatrix &R, 
                                         int num_clusters, Clusters &clusters,                                         
                                         const BrdfPointLights &all_vpls,
                                         const Pixels &pixels) {
    
    clusters.clear();
        
    if (all_vpls.size() < num_clusters) {
        Log::warn() << "Too few VPLs. Assign each VPL to a cluster." << endn;

        for (int i = 0; i < all_vpls.size(); ++i) {            
            Cluster new_cluster;
            new_cluster.center = i;    
            new_cluster.weight = 1.0f;
            clusters.push_back(new_cluster);
        }
        return;
    }

    Random &rd = *scene->get_random();                                         
    vector<Float> alphas;
    	
    // compute sum of distance cost from vertex i to all other vertices
    Stats stats;
    Log::info() << "Calculating distance..." << endn;
    stats.tic();
    R.generate_alphas_for_clustering(alphas);
    stats.toc();
    stats.print_elapsed_milliseconds();
    
    DiscretePdf pdf(alphas);

    Log::info() << "Clustering by sampling..." << endn;
    stats.tic();
    int num_cluster_attempts = 0;
    int count = 0;
    
    while (count < num_clusters) {
        num_cluster_attempts++;

        // sometimes it loops infinite when parameters are not properly set so here is for debug
        //if (num_cluster_attempts % 50 == 0)
        //    Log::info() << "Cluster attempts: " << num_cluster_attempts << "  count: " << count << endn;

        // this is the light index
        int k;
        Float pdf_k;
        pdf.sample(rd, k, pdf_k);
        
        if (pdf_k == 0.0f) continue;

        // only a few hundred or thousands of clusters, so just search
        int idx = -1;
        for (int i = 0; i < clusters.size(); ++i) {
            if (clusters[i].center == k) {
                idx = i;
                break;
            }
        }

        if (idx < 0) {
            Cluster new_cluster;
            new_cluster.center = k;    
            new_cluster.weight = 1.0f / pdf_k;
            clusters.push_back(new_cluster);
            count++;

        } else {
            
            Cluster &cluster = clusters[idx];
            cluster.weight += 1.0f / pdf_k;
        }
    }
    
    // take all selected columns into a separate matrix K
    // scale all selected columns by weight
    // note that the weight here is not for gathering. It is for clustering purpose only.
    LightTransportMatrix K(R.get_rows(), clusters.size());
    vector<Float> col(R.get_rows() * 3);
    for (int i = 0; i < clusters.size(); ++i) {
        int k = clusters[i].center;

        R.get_column(k, col);
        for (int j = 0; j < col.size(); ++j)
            col[j] *= clusters[i].weight;

        // this weight describes the reliability of the cluster center. When a column with low probability is chosen,
        // its weight is very high, and thus fewer columns (or might be none, even itself) are assigned to this cluster.

        K.set_column(i, col);
    }

    
    vector<Float> K_norms(K.get_cols());
    K.get_column_norms(K_norms);

    stats.toc();
    Log::info() << "Clustering by sampling ... [DONE]" << endn;
    stats.print_elapsed_milliseconds();
        
    // TODO:
    // If the min error is too high, we fork a new cluster. (This is a second step of exploration)
    // which does not require too much memory to build the clusters. 

    // It is possible to use the same method as computing alpha to compute the distance matrix
    // between all columns in R and all columns in K:
    // n_K^t * n_R - K^t * R
    //
    // NOTE: we do not implement the distance calculation as matrix multiplication here
    // because we might work with the reduced matrix R (which only uses a subset of all VPLs to compute matrix R)
    // this allows the clustering to scale better by avoiding establishing the full matrix R in memory

    // assign remaining columns to the nearest
    Log::info() << "Assigning columns ..." << endn;
    stats.tic();
    vector<Float> vpl_mags(all_vpls.size());

    vector<Float> mi(3 * pixels.size());
    vector<Float> mk(3 * pixels.size());
    for (int i = 0; i < all_vpls.size(); ++i) {
        // find the nearest cluster for VPL i
        int min_c = -1;
        Float min_dist = FLT_MAX;

        // this is light index
        R.get_column(i, mi);
        Float norm_mi = magnitude(mi);

        for (int j = 0; j < clusters.size(); ++j) {
            K.get_column(j, mk);
                        
            Float dist = fabs(K_norms[j] * norm_mi - dot(mi, mk));
            if (dist < min_dist) {
                min_c = j;
                min_dist = dist;
            }            
        }

        // assign column i to the nearest cluster
        clusters[min_c].indices.push_back(i);
        vpl_mags[i] = norm_mi;
    }
    stats.toc();
    Log::info() << "Assigning columns ... [DONE]" << endn;
    stats.print_elapsed_milliseconds();

    // generate probability for each columns in the cluster
    for (int i = 0; i < clusters.size(); ++i) {
        vector<Float> mags(clusters[i].indices.size());
        for (int j = 0; j < clusters[i].indices.size(); ++j) {
            int k = clusters[i].indices[j];

            mags[j] = vpl_mags[k];
        }
        clusters[i].pdf.set_distribution(mags);
    }

    // choose representative
    for (int i = 0; i < clusters.size(); ++i) {
        Cluster &c = clusters[i];

        int j;
        Float pdf_j;        
        clusters[i].pdf.sample(rd, j, pdf_j);
        
        if (j < 0) {
            if (c.indices.size() > 0)
                c.representative = c.indices[c.indices.size() / 2];
            else {
                c.representative = -1;
                // this case happens when the cluster center weight is too high.
            }

            c.weight = 0.0f;
        } else {        
            c.representative = c.indices[j];
            c.weight = 1.0f / pdf_j;
        }
    }

    // remove clusters with no representative
    Clusters old_clusters;
    old_clusters.assign(clusters.begin(), clusters.end());
    clusters.clear();
    for (int i = 0; i < old_clusters.size(); ++i) {
        if (old_clusters[i].representative >= 0) {
            clusters.push_back(old_clusters[i]);
        }
    }
}

} // end namespace
