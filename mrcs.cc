#include "mrcs.h"
#include "scene.h"
#include "camera.h"
#include "aggregate.h"
#include "shape.h"
#include "brdf_point_light.h"
#include "light_transport_matrix.h"
#include "cluster.h"
#include "cluster_builder.h"
#include "stats.h"
#include "vpl.h"
#include "log.h"

namespace Renzoku {

MatrixRowColumnSampling::MatrixRowColumnSampling() {
    max_clusters = 1000;
    block_size = 32;
    cluster_ratio = 0.01f;
}

/**
 * Hasan's SG07 row-column sampling approach
 * 
 * Gathering from cluster representatives.
 */
Rgb MatrixRowColumnSampling::gather_clusters(const Receiver &r, IVirtualPointLightEvaluator *vpl) {    
    if (r.m == NULL) return DefaultRgb::black;

    Rgb radiance;
    for (int i = 0; i < all_clusters.size(); ++i) {
        int k = all_clusters[i].representative;        
        if (k < 0) continue;

        BrdfPointLight &pl = (*all_vpls)[k];            
        Rgb Lo = vpl->radiance(scene, r, pl);        
        radiance += all_clusters[i].weight * Lo;        
	}
    return radiance;
}

void MatrixRowColumnSampling::initialize(Scene *scene, 
                                         IVirtualPointLightEvaluator *evaluator, 
                                         BrdfPointLights &vpls) {
    this->scene = scene;
    all_vpls = &vpls;
    num_clusters = std::max(1, (int)(cluster_ratio * all_vpls->size()));
    if (num_clusters > max_clusters) {
        num_clusters = max_clusters;
        Log::info() << "Too many clusters. Setting number of clusters to " << max_clusters << endn;
    }

    Random &rd = *scene->get_random();

    // generate pixels for reduced matrix
    Stats stats;
    Log::info() << "Generating image samples for reduced matrix..." << endn;
    stats.tic();
        ImageSampler img_sampler;
        Pixels pixels;
        img_sampler.sample_stratified(scene, block_size, pixels);
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "MRCS: " << pixels.size() << " x " << num_clusters << endn;

    // assume the number of clusters does not scale with number of VPLs, 
    // so have a cap at max VPLs that we should use for finding clustering.
    int max_columns = 10000;
    BrdfPointLights reduced_vpls;
    if (all_vpls->size() > max_columns) {
        vector<Float> power(all_vpls->size());
        for (int i = 0; i < all_vpls->size(); ++i) {
            power[i] = (*all_vpls)[i].power().avg();
        }
        DiscretePdf pdf(power);
                
        reduced_vpls.resize(max_columns);
        for (int i = 0; i < max_columns; ++i) {
            int k;
            Float pdf_k;
            pdf.sample(rd, k, pdf_k);
            reduced_vpls[i] = (*all_vpls)[k];
            reduced_vpls[i].set_contribution((*all_vpls)[k].power() / pdf_k / max_columns);
        }        
    } else {
        reduced_vpls.resize(all_vpls->size());
        for (int i = 0; i < all_vpls->size(); ++i) 
            reduced_vpls[i] = (*all_vpls)[i];        
    }

    // generate reduced matrix
    Log::info() << "Generating reduced matrix: " << pixels.size() << " x " << reduced_vpls.size() << endn;
    stats.tic();
        LightTransportMatrix R(pixels.size(), reduced_vpls.size());
        R.fill(scene, evaluator, reduced_vpls, pixels);
    stats.toc();
    stats.print_elapsed_milliseconds();

    /*
    cout << "Generating form factor matrix..." << endl;
    stats.tic();    
        FormFactorMatrix G;
        G.fill(scene, evaluator, *all_vpls, pixels);
    stats.toc();
    cout << "Elapsed: " << stats.elapsed() << endl;

    cout << "Generating visibility matrix..." << endl;
    stats.tic();    
        VisibilityMatrix V;
        V.fill(scene, evaluator, *all_vpls, pixels);
    stats.toc();
    cout << "Elapsed: " << stats.elapsed() << endl;

    cout << "Generating brdf light matrix..." << endl;
    stats.tic();    
        BrdfLightMatrix F1;
        F1.fill(scene, evaluator, *all_vpls, pixels);
    stats.toc();
    cout << "Elapsed: " << stats.elapsed() << endl;

    cout << "Generating brdf surface matrix..." << endl;
    stats.tic();    
        BrdfSurfaceMatrix F2;
        F2.fill(scene, evaluator, *all_vpls, pixels);
    stats.toc();
    cout << "Elapsed: " << stats.elapsed() << endl;

    ostringstream s1;
    s1 << scene->get_name() << "_mat_clamped.bin";
    R.save(s1.str().c_str());

    ostringstream s2;
    s2 << scene->get_name() << "_G.bin";
    G.save(s2.str().c_str());

    ostringstream s3;
    s3 << scene->get_name() << "_V.bin";
    V.save(s3.str().c_str());

    ostringstream s4;
    s4 << scene->get_name() << "_F1.bin";
    F1.save(s4.str().c_str());

    ostringstream s5;
    s5 << scene->get_name() << "_F2.bin";
    F2.save(s5.str().c_str());
    */

    // cluster columns of the reduced matrix
    Log::info() << "Clustering light transport matrix ..." << endn;    
    stats.tic();
        ClusterBuilder builder;    
        builder.cluster_by_sampling(scene, 
                                    evaluator,
                                    R, 
                                    num_clusters, all_clusters,                                     
                                    *all_vpls,
                                    pixels);
    stats.toc();
    Log::info() << "Clustering light transport matrix ... [DONE]" << endn;    
    stats.print_elapsed_milliseconds();
    
    /*
    // save clusters
    ofstream of("clusters_R.txt");
    for (int i = 0; i < all_clusters.size(); ++i) {
        of << "Cluster: " << i + 1 << endl;
        for (int j = 0; j < all_clusters[i].indices.size(); ++j) {
            int k = all_clusters[i].indices[j];

            of << k << " " << R.get_magnitude_of_column(k) << endl;
        }
    }
    of.close();    

    // save indices
    of.open("index_R.txt");
    for (int i = 0; i < all_clusters.size(); ++i) {
        for (int j = 0; j < all_clusters[i].indices.size(); ++j) {
            of << all_clusters[i].indices[j] << endl;
        }        
    }
    of.close();
    */
}

} // end namespace