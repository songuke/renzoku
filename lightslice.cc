#include "lightslice.h"
#include "scene.h"
#include "camera.h"
#include "aggregate.h"
#include "shape.h"
#include "concentric_mapping.h"
#include "brdf_point_light.h"
#include "cluster.h"
#include "cluster_builder.h"
#include "stats.h"
#include "vpl.h"
#include "log.h"

#include "image_block_view.h"

#include "mesh_view.h"  // debug only

#include <Eigen/Dense>
using namespace Eigen;

namespace Renzoku {

LightSlice::LightSlice() : kdtree(NULL), node_data(NULL), R(LightTransportMatrix(1, 1)) {

    num_clusters = 600;

    max_slice_size = 256;
    num_neighbor_slices = 16;
    max_slice_clusters = 2 * num_clusters;

    num_brdf_samples = 1024;
    
    incident_radiances.resize(num_brdf_samples);
    is_incident_radiance_reliable.resize(num_brdf_samples);
    brdf_samples.resize(num_brdf_samples);

    use_poly_interpolation = false;
    poly_deg = 5;

    // cone too small causes numerical error in the uniform sampling
    uniform_cone_half_angle = A_PI * 0.05f;

    // Metropolis sampling
    vector<Float> mutation_dist(Mutation::NUM_MUTATIONS);
    mutation_dist[Mutation::UNIFORM] = 0.33f;
    mutation_dist[Mutation::VPL] = 0.0f;    // this techinique, while biased towards VPLs, is not a problem as long as entire density over entire hemisphere is positive.
    mutation_dist[Mutation::CONE_PERTURB] = 0.33f;
    mutation_pdf.set_distribution(mutation_dist);

    density_radius = 0.1f;              // the radius controls bias
    num_sum_samples = 100000;

    perturb_cone_half_angle = A_PI * 0.05f;
    perturb_cos_theta_max = cos(perturb_cone_half_angle);
    perturb_cone_pdf = 1.0f / (TWO_PI * (1.0f - perturb_cos_theta_max));

    luminance_map_resolution = 32;      // map resolution in Jensen's method 

    same_local_wi = true;               // the sampled direction is the same in world or local space.
                                        // Apr 06:
                                        // same local: ensures direction always in upper hemisphere, and therefore avoid splotches
                                        // 
                                        // there is no real benefit between the two, but in terms of BRDF evaluation for opaque surface, 
                                        // is more efficient because at least the directions do not fail the under-surface test.

    view_similarity_cosine_max = cos(A_PI * 0.25f);  // this is a loose constraint. In most cases view difference, if not too severe, can still reuse the slice's map for sampling. 

    metropolis_neighbor_slices = 1;                 // by default for performance we just want to use the slice that contains the gather point
    metropolis_neighbor_slices_when_failed = 32;    // only when it's invalid then we look at its neighbors

    metropolis_tmp_slice_indices.reserve(metropolis_neighbor_slices_when_failed);
    metropolis_tmp_slice_pdf.reserve(metropolis_neighbor_slices_when_failed);    
    
    luminance_map_epsilon = 1e-3f;      // Jensen's approach: ratio of total illumination for each pixel, in which a direction is selected.
                                        
    density_epsilon_ratio = 1e-2f;      // relative to total incoming radiance by VPLs (without epsilon)

    distance_epsilon = 1.0f;           

    sampling_type = LocalSampler::METROPOLIS;
        
    radius_type = LocalDensity::FIXED_RADIUS;

    Log::info() << "LightSlice sampler type: " << sampling_type << endn;
}

bool LightSlice::is_view_similar(const Onb &uvn1, const Onb &uvn2) {    
      
    // FIXED: 
    // Apr 05
    //
    // View check causes splotches (patchy results) on curvy surfaces 
    //
    // A solution is to simply disable view check as we already cluster gather points into patches
    // A point that uses a sampling record that has view difference 
    // will still remains unbiased because all directions in the hemisphere
    // has a chance to be sampled.
    // Using such a record is still better than uniform sampling.
    //
    // Apr 06: 
    // We set the similarly angle to 45 deg, which is a very loose constraint. 
    // This avoids splotches, but also ensures more efficiency at extreme cases (where two patches face different directions)

    if (dot(uvn1.n(), uvn2.n()) < view_similarity_cosine_max) {
        return false;
    }

    if (same_local_wi) {
            
        if (dot(uvn1.u(), uvn2.u()) < view_similarity_cosine_max) {
            return false;
        }

    } else {
        
        // When we assume world wi is the same, it means we assume the tangent is the same, which leads to 
        // the local basis are already the same. 
        // Therefore, it is not necessary to care about right vector at the receiver.
        //
        // But this will not support anisotropic materials well unless quads are used.

    }
        
    return true;
}

inline static Vec6 concat_vec6(const Vec3 &u, const Vec3 &v) {
    Vec6 p;    
    p.p[0] = u.x();
    p.p[1] = u.y();
    p.p[2] = u.z();
    p.p[3] = v.x();
    p.p[4] = v.y();
    p.p[5] = v.z();
    return p;
}

/**
 * Mapping from vertex position and normal to 6D space
 * as suggested in Multidimensional Lightcuts paper.
 */
inline static Vec6 map_vec6(const Vec3 &v, const Vec3 &n, Material *m, const Vec3 &wo, Float diagonal) {
    
    const Float weight = 1.0f / 16.0f * diagonal;
    const Float max_gloss_scale = 4.0f;
    const Float max_gloss_exponential = 1000.0f;
        
    Vec3 directionality;

    if (m) {
        switch (m->get_bsdf()->get_bsdf_type()) {

            case Bsdf::PHONG:
            case Bsdf::WARD:
            {
                ISpecular *bsdf = dynamic_cast<ISpecular *>(m->get_bsdf());
                Vec3 wr = 2.0f * dot(wo, n) * n - wo;
                Float gloss_scale = std::min(max_gloss_scale, std::max(1.0f, bsdf->get_gloss_exponential() / max_gloss_exponential));
                directionality = wr * gloss_scale;
                break;
            }

            case Bsdf::MIRROR:
            case Bsdf::GLASS:
            {
                Vec3 wr = 2.0f * dot(wo, n) * n - wo;        
                directionality = wr * max_gloss_scale;
                break;
            }

            case Bsdf::LAMBERTIAN:
            default:
                directionality = n;
                break;
        }
    } else {
        directionality = n;
    }

    Vec3 position = v;

    Vec6 p;
    p.p[0] = position.x();
    p.p[1] = position.y();
    p.p[2] = position.z();
    p.p[3] = directionality.x() * weight;
    p.p[4] = directionality.y() * weight;
    p.p[5] = directionality.z() * weight;
    return p;
}

inline static Vec9 map_vec9(const Vec3 &v, const Vec3 &n, const Vec3 &u, Material *m, const Vec3 &wo, Float diagonal) {
    // lower weight ensures smaller patch but normal and tangent constraint becomes loose.
    // increase the number of slices if we want smaller patch but keep normal and tangent constraint.
    const Float weight = 1.0f / 16.0f * diagonal;
    const Float max_gloss_scale = 4.0f;
    const Float max_gloss_exponential = 1000.0f;
        
    Vec3 directionality;
    Vec3 horizontal;

    if (m) {
        switch (m->get_bsdf()->get_bsdf_type()) {

            case Bsdf::PHONG:
            case Bsdf::WARD:
            {
                ISpecular *bsdf = dynamic_cast<ISpecular *>(m->get_bsdf());
                Vec3 wr = 2.0f * dot(wo, n) * n - wo;
                Float gloss_scale = std::min(max_gloss_scale, std::max(1.0f, bsdf->get_gloss_exponential() / max_gloss_exponential));
                directionality = wr * gloss_scale;

                horizontal = u * gloss_scale;
                break;
            }

            case Bsdf::MIRROR:
            case Bsdf::GLASS:
            {
                Vec3 wr = 2.0f * dot(wo, n) * n - wo;        
                directionality = wr * max_gloss_scale;
                horizontal = u * max_gloss_scale;
                break;
            }

            case Bsdf::LAMBERTIAN:
            default:
                directionality = n;
                horizontal = u;
                break;
        }
    } else {
        directionality = n;
        horizontal = u;
    }

    Vec3 position = v;

    Vec9 p;
    p.p[0] = position.x();
    p.p[1] = position.y();
    p.p[2] = position.z();
    p.p[3] = (directionality.x()) * weight;
    p.p[4] = (directionality.y()) * weight;
    p.p[5] = (directionality.z()) * weight;
    p.p[6] = horizontal.x() * weight;
    p.p[7] = horizontal.y() * weight;
    p.p[8] = horizontal.z() * weight;
    return p;
}

void LightSlice::generate_slices() {
    Size2 img_size = scene->get_image_size();
    int num_pixels = img_size.height * img_size.width;

    if (node_data) delete [] node_data;
    node_data = new NodeData6[num_pixels];
    int num_nodes = 0;

    // create shading points    
    Camera *camera = scene->get_camera();
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float tick = scene->get_tick();
    Float max_tmax = scene->get_max_tmax();    
    for (int i = 0; i < img_size.height; ++i) {
        for (int j = 0; j < img_size.width; ++j) {
            Ray r = camera->shoot_ray(i, j);

            HitRecord rec;
            // make sure material is available
            if (agg->hit(r, tmin, max_tmax, tick, rec) && rec.material) {
                Receiver recv(r, rec);                
                node_data[num_nodes].p = map_vec6(recv.p, recv.shading_n, recv.m, recv.wo, scene_bb_diag);
                node_data[num_nodes].index = i * img_size.width + j;
                num_nodes++;
            }
        }
    }
    
    Log::info() << "Surface gather points: " << num_nodes << endn;
    
    // TODO: generate more surface gather points (not from camera) to support multi-bounce tracing

    Log::info() << "Building kd-tree for gather points..." << endn;
    vector<pair<int, int>> slice_clusters;
    if (kdtree) delete kdtree;
    kdtree = new KdTree6(max_slice_size, node_data, num_nodes, scene, slice_clusters);

    Log::info() << "Slices: " << slice_clusters.size() << endn;
    
    Random &rd = *scene->get_random();
    slices.clear();
    slices.resize(slice_clusters.size());
    for (int i = 0; i < slice_clusters.size(); ++i) {
        slices[i].start = slice_clusters[i].first;
        slices[i].end = slice_clusters[i].second;
                        
        slices[i].representative = slices[i].start + rd() * (slices[i].end - slices[i].start + 1);        

        slices[i].color = Rgb::from_hsv(rd() * 360.0f, 0.75f, 1.0f);    
    }
}


static void minus(const vector<Float> &a, const vector<Float> &b, vector<Float> &ab) {
    ab.resize(a.size());
    for (int i = 0; i < a.size(); ++i) 
        ab[i] = a[i] - b[i];
}

static Float dot(const vector<Float> &a, const vector<Float> &b) {
    Float d = 0.0f;
    for (int i = 0; i < a.size(); ++i)
        d += a[i] * b[i];
    return d;
}

static Float sum(const vector<Float> &a) {
    Float s = 0.0f;
    for (int i = 0; i < a.size(); ++i)
        s += a[i];
    return s;
}

static void pick(const vector<Float> &in, const vector<int> &indices, vector<Float> &out) {
    out.resize(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        out[i] = in[indices[i]];
    }
}

/**
 * Compute distance from each column to every column in the matrix
 */
static Float compute_distance(MatrixC<Float> &R, vector<Float> &column_norms,
                              MatrixC<Float> &Rt,
                              vector<Float> &temp_col, vector<Float> &temp_row) {
    int rows = R.get_rows();
    int cols = R.get_cols();
    
    // alpha = (n_R^t * n_R - R^t * R) 1
    //       = (n_R^t * n_R * 1 - R^t * R * 1)
            
    // temp_col u has rows elements
    // temp_row v has cols elements
    R.row_sum(temp_col);                       // u = R * 1    
    Rt.mul(temp_col, temp_row);                // v = Rt * u
    
    Float c = sum(column_norms);        // n_R * 1
    
    Float cost = 0.0f;
    for (int i = 0; i < cols; ++i) {
        cost += column_norms[i] * c - temp_row[i];
    }
    return cost;
}

static Float compute_distance_pos(Cluster &cluster, const BrdfPointLights &vpls) {
    Float cost = 0.0f;
    for (int i = 0; i < cluster.indices.size(); ++i) {
        for (int j = 0; j < cluster.indices.size(); ++j) {
            cost += (vpls[cluster.indices[i]].org() - vpls[cluster.indices[j]].org()).squared_length();
        }
    }
    return cost;
}


struct ClusterCost {
    int index;
    Float cost;

    ClusterCost(int index, Float cost) : index(index), cost(cost) {}
};

struct ClusterCostLess {    
    bool operator()(ClusterCost &a, ClusterCost &b) const {
        return a.cost < b.cost;
    }
};

typedef priority_queue< ClusterCost, vector<ClusterCost>, ClusterCostLess >  ClusterCostHeap;

static void choose_representative(Random &rd, Cluster &c) {
    int j;
    Float pdf_j;
    c.pdf.sample(rd, j, pdf_j);

    if (j < 0) {
        // NOTE: all norms are zero. Just use the median.
        c.representative = c.indices[c.indices.size() / 2];
        c.weight = 0.0f;
    }
    else {
        c.representative = c.indices[j];
        c.weight = 1.0f / pdf_j;
    }
}

static bool pair_less(const pair<int, Float> &a, const pair<int, Float> &b) {
    return a.second < b.second;
}

static void split_cluster(Random &rd, const Cluster &cluster, const MatrixC<Float> &matC,
                          vector<int> &indices1, vector<int> &indices2) {
                              
    indices1.clear();
    indices2.clear();
    if (matC.get_cols() == 1) {
        indices1.push_back(cluster.indices[0]);
        return;
    }

    if (matC.get_cols() == 2) {
        indices1.push_back(cluster.indices[0]);
        indices2.push_back(cluster.indices[1]);
        return;
    }
    
    vector<Float> a, b;
    vector<Float> minus_ab;
    int i, j;
    Float pdf_i, pdf_j;

    cluster.pdf.sample(rd, i, pdf_i);
    cluster.pdf.sample(rd, j, pdf_j);

    if (i < 0 || j < 0) {
        Log::info() << "Split by half due to sampling distribution problems." << endn;

        vector<int>::const_iterator mid = cluster.indices.begin() + cluster.indices.size() / 2;
        indices1.assign(cluster.indices.begin(), mid);
        indices2.assign(mid, cluster.indices.end());

        return;
    }

    matC.get_column(i, a);
    matC.get_column(j, b);
    minus(a, b, minus_ab);

    vector<Float> c;
    vector<Float> minus_ac;

    vector<pair<int, Float>> line;
    line.reserve(matC.get_cols());
    for (int k = 0; k < matC.get_cols(); ++k) {
        if (k == i || k == j) continue;

        matC.get_column(k, c);

        minus(a, c, minus_ac);
                                                // Float t = dot(minus_ac, minus_ab) / sqr_ab;
        Float t = dot(minus_ac, minus_ab);      // no need to divide by square of ab since it does not vary among points

        line.push_back(pair<int, Float>(k, t));
    }

    std::sort(line.begin(), line.end(), pair_less);
    
    Float max_dist = 0.0f;
    int max_k = -1;
    for (int k = 1; k < line.size(); ++k) {
        Float dist = line[k].second - line[k - 1].second;
        if (dist > max_dist) {
            max_dist = dist;
            max_k = k;
        }
    }

    if (max_k == -1) {        
        
        // just split by half and return because of degeneracy: ab is not a line, or the dot product is always zero        
        vector<int>::const_iterator mid = cluster.indices.begin() + cluster.indices.size() / 2;
        indices1.assign(cluster.indices.begin(), mid);
        indices2.assign(mid, cluster.indices.end());
        
    } else {

        indices1.push_back(cluster.indices[i]);
        indices2.push_back(cluster.indices[j]);
        for (int k = 0; k < max_k; ++k)
            indices1.push_back(cluster.indices[line[k].first]);
        for (int k = max_k; k < line.size(); ++k)       // max_k must be >= 0 otherwise compare -1 with unsigned int will result in false
            indices2.push_back(cluster.indices[line[k].first]);
    }
}

static void split_cluster_pos(Random &rd, const Cluster &cluster, const BrdfPointLights &vpls,
                              Float scene_bb_diag,
                              vector<int> &indices1, vector<int> &indices2) {
                              
    indices1.clear();
    indices2.clear();
    if (cluster.indices.size() == 1) {
        indices1.push_back(cluster.indices[0]);
        return;
    }

    if (cluster.indices.size() == 2) {
        indices1.push_back(cluster.indices[0]);
        indices2.push_back(cluster.indices[1]);
        return;
    }
    
    // map everything into 6D, find the longest axis, and split into two.
    const int dim = 6;
    int size = cluster.indices.size();

    vector<NodeData6> points(size);
    for (int i = 0; i < size; ++i) {
        const BrdfPointLight &vpl = vpls[cluster.indices[i]];
        points[i].p = map_vec6(vpl.org(), vpl.normal(), vpl.get_material(), vpl.get_wi(), scene_bb_diag);
        points[i].index = cluster.indices[i];
    }
    
    // take node's bounding box and determine the longest axis
    BoundingBoxN<dim> bb;
    for (int i = 0; i < size; ++i)
        bb.merge(points[i].p);
    int axis = bb.get_longest_axis();

    // sort the points on the axis
    struct ComparePoint {
        int axis;

        ComparePoint(int axis) : axis(axis) {}

        bool operator()(const NodeData<dim> &a, const NodeData<dim> &b) {
            return a.p[axis] < b.p[axis];
        }
    };
    ComparePoint compare_point(axis);
    std::sort(points.begin(), points.end(), compare_point);

    // and split at the longest segment    
    Float max_dist = 0.0f;
    int max_i = -1;
    for (int i = 1; i < size; ++i) {
        Float dist = points[i].p[axis] - points[i - 1].p[axis];        
        if (dist > max_dist) {
            max_dist = dist;
            max_i = i;
        }
    }
    
    if (max_i == -1) {        
        
        // just split by half and return 
        vector<int>::const_iterator mid = cluster.indices.begin() + cluster.indices.size() / 2;
        indices1.assign(cluster.indices.begin(), mid);
        indices2.assign(mid, cluster.indices.end());
        
    } else {

        for (int i = 0; i < max_i; ++i)
            indices1.push_back(points[i].index);
        for (int i = max_i; i < size; ++i)
            indices2.push_back(points[i].index);
    }
}

static FILE *open_ply(const string &name, int num_vertices) {
    FILE *f = fopen(name.c_str(), "w");
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "element vertex %d\n", num_vertices);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "property uchar red\n");
    fprintf(f, "property uchar green\n");
    fprintf(f, "property uchar blue\n");
    fprintf(f, "end_header\n");    
    return f;
}

static void close_ply(FILE *f) {
    fclose(f);
}

static void write_ply(FILE *f, const Vec3 &p, const Rgb &color) {    
    fprintf(f, "%f %f %f %d %d %d\n", 
        p.x(), p.y(), p.z(), 
        (int)(color.red() * 255), (int)(color.green() * 255), (int)color.blue() * 255);    
}

void LightSlice::initialize(Scene *scene, 
                            IVirtualPointLightEvaluator *evaluator, 
                            BrdfPointLights &vpls) {
        
    Stats init_stats;
    init_stats.tic();

    this->scene = scene;
    this->scene_bb_diag = scene->get_bounding_box().diagonal();
    this->evaluator = evaluator;

    if (scene->get_image_view()) {
        scene->get_image_view()->mouse_event.attach_observer(this);
    }

    all_vpls = &vpls;

    // HACK: currently we only work on 2-bounce path tracing due to the slice cluster. Therefore, we should only work with VPL that has bounce at most once.
    // guiding directional tracing with deeper bounce VPL should only be used when we are not confined to only able to sample at second bounce.
    BrdfPointLights one_bounce_vpls;
    for (int i = 0; i < vpls.size(); ++i)
        if (vpls[i].get_bounce() <= 1)
            one_bounce_vpls.push_back(vpls[i]);
    vpls.clear();
    vpls.assign(one_bounce_vpls.begin(), one_bounce_vpls.end());
    // this hack should be removed when our cache points can be used at multiple bounces


    if (vpls.size() == 0) {
        Log::info() << "There are no VPLs to cluster." << endn;
        return;
    }
    
    Random &rd = *scene->get_random();
    
    // slicing
    Stats stats;
    Log::info() << "Generating slices ..." << endn;
    stats.tic();

    generate_slices();
    
    /*
    MeshView *mesh_view = scene->get_mesh_view();
    BoundingBoxes bbs;
    kdtree->get_cluster_bounding_boxes(bbs);
    mesh_view->set_bounding_boxes(bbs);
    return;
    */

    sampling_records.resize(slices.size());
    
    // test nearest clusters
    /*
    Vec9 p = concat_vec6(Vec3(-277.762299, 272.726135, -559.417358), Vec3(0.000000000, 0.191311017, 240.248581));

    vector<int> slice_indices;
    kdtree->find_nearest_cluster(p, 16, slice_indices);    
    
    for (int i = 0; i < slice_indices.size(); ++i) {
        int slice_index = slice_indices[i];
        slices[slice_index].color = DefaultRgb::red * (Float)(i+1) / slice_indices.size();
    }*/
            
    Pixels pixels(slices.size());
    Size2 img_size = scene->get_image_size();
    for (int i = 0; i < slices.size(); ++i) {
        int selected = slices[i].representative;
        
        int pixel_index = node_data[selected].index;
        int y = pixel_index / (int)img_size.width;
        int x = pixel_index % (int)img_size.width;

        pixels[i] = Vec2(x, y);
    }
        
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "Tentative reduced matrix: " << pixels.size() << " x " << num_clusters << endn;

    // assume the number of clusters does not scale with number of VPLs, 
    // so have a cap at max VPLs that we should use for finding clustering.    
    /*
    int max_columns = 1000000;
    BrdfPointLights reduced_vpls;
    if (all_vpls->size() > max_columns) {
        Log::info() << "Too large reduced matrix. Downsampling to " << max_columns << " columns." << endn;
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
            reduced_vpls[i].set_throughput((*all_vpls)[k].power() / pdf_k / max_columns);
        }        
    } else {
        reduced_vpls.resize(all_vpls->size());
        for (int i = 0; i < all_vpls->size(); ++i) 
            reduced_vpls[i] = (*all_vpls)[i];        
    }
    */

    // generate reduced matrix
    Log::info() << "Generating reduced matrix: " << pixels.size() << " x " << vpls.size() << endn;
    stats.tic();

    R.resize(pixels.size(), vpls.size());    
    MatrixC<Rgb> &matR = R.get_matrix();
    
    //matR.resize(pixels.size(), vpls.size());
    //matI.resize(pixels.size(), vpls.size());
        
    matV.resize(pixels.size(), vpls.size());
    
    // initialize matrix entries
    matR.set(0.0f);
    matV.set(0.0f);
    
    Camera *camera = scene->get_camera();
    Aggregate *agg = scene->get_aggregate();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Lights &lights = *scene->get_lights();
    
    for (int ridx = 0; ridx < pixels.size(); ++ridx) {
        Vec2 pixel = pixels[ridx];
        
        Ray r = camera->shoot_ray(pixel.y(), pixel.x());
        HitRecord rec;
        if (! agg->hit(r, tmin, max_tmax, tick, rec) || rec.material == NULL) {
            Log::info() << "Bad slice representative: " << ridx << endn;
            continue;
        }

        Receiver recv(r, rec);
        for (int k = 0; k < vpls.size(); ++k) {                        
            bool visible = (agg->hit(recv.p, vpls[k].org(), tmin, tick) == false);
            //Rgb radiance = evaluator->radiance(scene, recv, vpls[k], visible);

            // NOTE: should use radiance matrix instead of incident radiance matrix because
            // + The matrix clustering is used for out-going radiance evaluation. Clustering with incident radiance causes artifacts in output. 
            
            // only cluster incident radiance
            Rgb radiance = (Float)visible * evaluator->incident_radiance(scene, recv, vpls[k]);
            
            matR(ridx, k) = radiance;                 
            matV(ridx, k) = (Float)visible;            
        }        

        // keep uvn for Metropolis sampling record
        sampling_records[ridx].uvn = recv.shading_uvn;
    }

    stats.toc();
    stats.print_elapsed_milliseconds();

    // cluster columns of the reduced matrix
    Log::info() << "Clustering light transport matrix..." << endn;    
    stats.tic();
        ClusterBuilder builder;    
        builder.cluster_by_sampling(scene, 
                                    evaluator,
                                    R, 
                                    num_clusters, mrcs_clusters,
                                    *all_vpls,
                                    pixels);
    stats.toc();
    stats.print_elapsed_milliseconds();
    Log::info() << "MRCS clusters: " << mrcs_clusters.size() << endn;
    

    /*
    // NOTE: only for testing. Fake MRCS clusters by making each cluster a VPL
    Log::info() << "Using a fake MRCS clustering..." << endn;
    max_slice_clusters = all_vpls->size();
    vector<float> fake_pdf;
    fake_pdf.push_back(1.0f);
    mrcs_clusters.resize(all_vpls->size());
    for (int i = 0; i < all_vpls->size(); ++i) {
        mrcs_clusters[i].weight = 1.0f;
        mrcs_clusters[i].indices.clear();
        mrcs_clusters[i].indices.push_back(i);
        mrcs_clusters[i].representative = i;
        mrcs_clusters[i].pdf.set_distribution(fake_pdf);
    }*/

    // slice refinement
    Log::info() << "Slice clustering..." << endn;

    Log::info() << "Slice memory allocating..." << endn;
    stats.tic();

    int max_cols = 0;
    for (int i = 0; i < mrcs_clusters.size(); ++i) {
        if (max_cols < mrcs_clusters[i].indices.size())
            max_cols = mrcs_clusters[i].indices.size();
    }

    // pre-allocate some temporary memory
    vector<int> slice_indices;
    slice_indices.reserve(num_neighbor_slices);     

    MatrixC<Rgb> matL;                              // matrix that extracts rows from R
    matL.reserve(num_neighbor_slices, max_cols);

    MatrixC<Rgb> matC;                              // matrix that extracts rows and cols from R
    matC.reserve(num_neighbor_slices, max_cols);
    
    MatrixC<Float> matCf;                           // flatten version of L
    matCf.reserve(num_neighbor_slices * 3, max_cols);
        
    MatrixC<Float> matCfTrans;
    matCfTrans.reserve(max_cols, num_neighbor_slices * 3);
    
    matrix_clusters.clear();
    matrix_clusters.resize(slices.size());          // store clusters per slice
    for (int i = 0; i < slices.size(); ++i) {       // memory preallocation - indices still might grow
        matrix_clusters[i].reserve(max_slice_clusters);
    }

    vector<Float> norms;    
    vector<int> indices1;
    vector<int> indices2;

    vector<Float> temp_col;
    //temp_col.resize(num_neighbor_slices * 3);
    vector<Float> temp_row;
    //temp_row.resize(max_cols);

    vector<Float> norms_R;
    R.get_column_norms(norms_R);

    stats.toc();
    stats.print_elapsed_milliseconds();
        
    Log::info() << "Refining clustering per slice..." << endn;
    stats.tic();

    global_incident_cache.clear();
    
    for (int s = 0; s < slices.size(); ++s) {
        
        // find neighbor slices and extract the slice matrix L

        // take a first point in the slice and look for neighbors
        Vec6 &p = node_data[slices[s].representative].p;        
        kdtree->find_nearest_cluster(p, num_neighbor_slices, slice_indices);    
    
        // extract rows from matrix R        
        matR.subrow(slice_indices, matL);
        
        // determine cost for each current cluster by L and push to heap        
        Clusters &slice_clusters = matrix_clusters[s];
        slice_clusters.assign(mrcs_clusters.begin(), mrcs_clusters.end());

        
        ClusterCostHeap heap;        
        for (int c = 0; c < slice_clusters.size(); ++c) {

            // take all columns that belong to the current cluster                        
            matL.subcol(slice_clusters[c].indices, matC);
            matC.flatten(matCf);
            matCf.get_column_norm(norms);
            matCf.transpose(matCfTrans);

            Float cost = compute_distance(matCf, norms, matCfTrans, temp_col, temp_row);
            heap.push(ClusterCost(c, cost));

            // TODO: test if the cost is the same as brute-force            
        }
                
        
        while (heap.size() > 0 && slice_clusters.size() < max_slice_clusters) {
            // take the most expensive cluster
            ClusterCost cc = heap.top();
            
            // and split it into two
            if (slice_clusters[cc.index].indices.size() <= 1) break;

            matL.subcol(slice_clusters[cc.index].indices, matC);
            matC.flatten(matCf);
            split_cluster(rd, slice_clusters[cc.index], matCf, indices1, indices2);            

            // make two new clusters
            Cluster &cluster1 = slice_clusters[cc.index];
            
            matL.subcol(indices1, matC);            
            matC.flatten(matCf);
            matCf.get_column_norm(norms);        // since we don't have the full in-memory R
            matCf.transpose(matCfTrans);
            Float cost1 = compute_distance(matCf, norms, matCfTrans, temp_col, temp_row);

            cluster1.indices.assign(indices1.begin(), indices1.end());
            pick(norms_R, indices1, norms);      // use norms of R for representative sampling
            cluster1.pdf.set_distribution(norms);
            choose_representative(rd, cluster1);


            Cluster cluster2;

            matL.subcol(indices2, matC);            
            matC.flatten(matCf);
            matCf.get_column_norm(norms);
            matCf.transpose(matCfTrans);
            Float cost2 = compute_distance(matCf, norms, matCfTrans, temp_col, temp_row);

            cluster2.indices.assign(indices2.begin(), indices2.end());
            pick(norms_R, indices2, norms);
            cluster2.pdf.set_distribution(norms);
            choose_representative(rd, cluster2);

            // and insert new costs to heap
            int new_index = slice_clusters.size();
            ClusterCost cc1(cc.index, cost1);
            ClusterCost cc2(new_index, cost2);
            heap.pop();
            heap.push(cc1);
            heap.push(cc2);

            slice_clusters.push_back(cluster2);
        }
        
        // our method of cluster refinement based on positions
        /*
        for (int c = 0; c < slice_clusters.size(); ++c) {
            Float cost = compute_distance_pos(slice_clusters[c], vpls);
            heap.push(ClusterCost(c, cost));
        }
                
        while (heap.size() > 0 && slice_clusters.size() < max_slice_clusters) {
            // take the most expensive cluster
            ClusterCost cc = heap.top();
            
            // and split it into two
            if (slice_clusters[cc.index].indices.size() <= 1) break;
            split_cluster_pos(rd, slice_clusters[cc.index], vpls, scene_bb_diag, indices1, indices2);

            // make two new clusters
            Cluster &cluster1 = slice_clusters[cc.index];
            cluster1.indices.assign(indices1.begin(), indices1.end());
            pick(norms_R, indices1, norms);      // use norms of R for representative sampling
            cluster1.pdf.set_distribution(norms);
            choose_representative(rd, cluster1);
            Float cost1 = compute_distance_pos(cluster1, vpls);

            Cluster cluster2;
            cluster2.indices.assign(indices2.begin(), indices2.end());
            pick(norms_R, indices2, norms);
            cluster2.pdf.set_distribution(norms);
            choose_representative(rd, cluster2);
            Float cost2 = compute_distance_pos(cluster2, vpls);

            // and insert new costs to heap
            int new_index = slice_clusters.size();
            ClusterCost cc1(cc.index, cost1);
            ClusterCost cc2(new_index, cost2);
            heap.pop();
            heap.push(cc1);
            heap.push(cc2);

            slice_clusters.push_back(cluster2);
        }*/


        // compute average incident radiance and cone bound for each cluster
        
        // get receiver
        // IMPROVE: in fact receiver is not needed. Only p and n is needed to compute the incident radiance so this can be cached in node data.
        int pixel_index = node_data[slices[s].representative].index;
        int y = pixel_index / (int)img_size.width;
        int x = pixel_index % (int)img_size.width;
        Ray ray = camera->shoot_ray(y, x);
        HitRecord rec;
        if (! agg->hit(ray, tmin, max_tmax, tick, rec) || rec.material == NULL) {
            Log::info() << "Error. Ray should hit and material not null." << endn;
        }        

        Receiver recv(ray, rec);

        vector<IncidentLight> cache;

        for (int c = 0; c < slice_clusters.size(); ++c) {
            
            Cluster &cluster = slice_clusters[c];

            Rgb incident;
            Rgb outgoing;
            for (int i = 0; i < cluster.indices.size(); ++i) {
                BrdfPointLight &vpl = vpls[cluster.indices[i]];
                
                Float visibility = matV(s, cluster.indices[i]);
                Rgb incident_power = evaluator->incident_radiance(scene, recv, vpl);
                
                incident += visibility * incident_power;                
                outgoing += matR(s, cluster.indices[i]);

                // only support opaque surface at this moment
                Vec3 wi = unit_vector(vpl.org() - recv.p);
                if (visibility > 0.0f && dot(recv.shading_n, wi) > ZERO_EPSILON && incident_power != DefaultRgb::black) {
                    cluster.bc.merge_direction(wi);
                    //cluster.dirs.push_back(wi);

                }

                // build global incident cache
                if (cluster.indices[i] == cluster.representative && visibility > 0.0f) {
                    cache.push_back(IncidentLight(wi, incident_power, c));
                }

                // store wi
                if (cluster.indices[i] == cluster.representative) {
                    cluster.slice_wi = wi;

                    // upper hemisphere means visible. We need to check the power because of the BRDF at the VPL.
                    //
                    // we don't care about r.wo, so do not use pdf(recv, wi) function which might check if r.wo is under surface
                    if (visibility > 0.0f && dot(recv.shading_n, wi) > ZERO_EPSILON && incident_power != DefaultRgb::black) {                    
                        cluster.is_upper_hemisphere = true;
                    }
                }
            }

            // CAN FIX:
            // just to be safe, make sure the bounding cone is capped at half angle HALF_PI
            if (cluster.bc.half_angle() > HALF_PI)
                cluster.bc.set_half_angle(HALF_PI);

            // The cone merge algorithm is not optimal.
            // Therefore, it might produce a reflex cone.
            //if (cluster.bc.half_angle() > HALF_PI) 
            /*
            {
                int a = 1;

                Cone tmp(Vec3(0.0f), 0.0f);
                    
                Cone tmp2(Vec3(0.0f), 0.0f);
                for (int k = 0; k < cluster.dirs.size(); ++k) {
                    if (dot(cluster.dirs[k], recv.shading_n) <= ZERO_EPSILON) {
                        Log::info() << "Dir below surface." << endn;
                    }

                    Log::info() << "dir: " << cluster.dirs[k] << endn;

                    if (k == 23) {
                        int a = 1;
                    }

                    tmp.merge_direction(cluster.dirs[k]);                        
                    Log::info() << "k = " << k << " Cone 1: " << tmp.normal() << "\t" << tmp.half_angle() << endn;

                    tmp2.merge_direction_algebraic(cluster.dirs[k]);
                    Log::info() << "k = " << k << " Cone 2: " << tmp2.normal() << "\t" << tmp2.half_angle() << endn;
                }
            }
            */

            cluster.total_incident_radiance = incident;
            cluster.total_outgoing_radiance = outgoing;
        }  

        
        global_incident_cache.push_back(cache);

        /*
        // NOTE: this is clustering for a particular slice!
        if (s == 215) {
            // write clusters into a PLY file with colors
            FILE *f = open_ply("representative.ply", slice_clusters.size());
            for (int c = 0; c < slice_clusters.size(); ++c) {

                Rgb color = Rgb::from_hsv(rd() * 360.0f, 0.75f, 1.0f);
                ostringstream oss;
                oss << "cluster" << c << ".ply";
                string cluster_name = oss.str();
                FILE *f2 = open_ply(cluster_name, slice_clusters[c].indices.size());
                for (int k = 0; k < slice_clusters[c].indices.size(); ++k) {
                
                    int idx = slice_clusters[c].indices[k];

                    Vec3 p = vpls[idx].org();

                    if (idx == slice_clusters[c].representative) {
                        write_ply(f, p, color);
                    }

                    write_ply(f2, p, color);
                }
                close_ply(f2);
            }
            close_ply(f);
        }
        */

        /*
        Float min_angle = A_PI;
        Float max_angle = 0.0f;
        for (int j = 0; j < slice_clusters.size(); ++j) {
            Cluster &cluster = slice_clusters[j];
            if (cluster.bc.half_angle() == 0.0f) {
                int a = 1;
            }
            
            if (cluster.bc.half_angle() > 0.0f && cluster.bc.half_angle() < min_angle)
                min_angle = cluster.bc.half_angle();
            if (cluster.bc.half_angle() > max_angle)
                max_angle = cluster.bc.half_angle();
        }
        Log::info() << "Min/max angle: " << min_angle << "\t" << max_angle << endn;
        */

        

    }    // end slice
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "Cache cluster norms for pdfs..." << endn;
    stats.tic();

    // cache cluster norm for each slice
    vector<int> reps(max_slice_clusters);
    for (int s = 0; s < slices.size(); ++s) {

        // take a first point in the slice and look for neighbors
        Vec6 &p = node_data[slices[s].representative].p;        
        kdtree->find_nearest_cluster(p, num_neighbor_slices, slice_indices);    
    
        // extract rows from matrix R        
        matR.subrow(slice_indices, matL);
        
        reps.resize(matrix_clusters[s].size());
        for (int c = 0; c < matrix_clusters[s].size(); ++c) {
            reps[c] = matrix_clusters[s][c].representative;
        }
        
        // take all representative columns that belong to the current cluster                        
        matL.subcol(reps, matC);

        vector<Float> cluster_norm;
        matC.get_column_norm(cluster_norm);

        // estimate cone size based on nearest direction (use kd-tree)
        /*
        Vec3 pos = Vec3(p[0], p[1], p[2]);
        for (int c = 0; c < matrix_clusters[s].size(); ++c) {
            BrdfPointLight &vpl1 = (*all_vpls)[matrix_clusters[s][c].representative];            
            Vec3 dir1 = unit_vector(vpl1.org() - pos);

            Float max_cosine = -1.0f;            
            for (int d = 0; d < matrix_clusters[s].size(); ++d) {
                BrdfPointLight &vpl2 = (*all_vpls)[matrix_clusters[s][d].representative];
                Vec3 dir2 = unit_vector(vpl2.org() - pos);

                Float cosine = dot(dir1, dir2);
                if (cosine > max_cosine) {
                    max_cosine = cosine;
                }
            }

            Float half_angle = acos(max_cosine);
            matrix_clusters[s][c].bc = Cone(dir1, half_angle);
        }*/

        // disable sampling if 
        // - is in lower hemisphere
        // - cone half angle is larger than PI/2
        for (int c = 0; c < matrix_clusters[s].size(); ++c) {

            // must take weight into account
            cluster_norm[c] *= matrix_clusters[s][c].weight;

            if (matrix_clusters[s][c].is_upper_hemisphere == false)
                cluster_norm[c] = 0.0f;

            //if (matrix_clusters[s][c].bc.half_angle() >= HALF_PI)
            //    cluster_norm[c] = 0.0f;
        }

        slices[s].cluster_pdf.set_distribution(cluster_norm);
    }
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "Make kd-tree for each slice..." << endn;
    stats.tic();

    // make kd-tree for each slice
    for (int s = 0; s < matrix_clusters.size(); ++s) {
        
        int num_visible_clusters = 0;
        
        int slice = s;
        for (int c = 0; c < matrix_clusters[slice].size(); ++c) {
            if (matrix_clusters[slice][c].is_upper_hemisphere) {
                // approximate and conservative upperbound
                num_visible_clusters++; 
            }
        }

        slices[s].tree_node_data = new NodeData2[num_visible_clusters];
        NodeData2 *tree_node_data = slices[s].tree_node_data;

        num_visible_clusters = 0;
                        
        if (! sampling_records[slice].is_valid) continue;

        for (int c = 0; c < matrix_clusters[slice].size(); ++c) {
            if (matrix_clusters[slice][c].is_upper_hemisphere) {
                                                            
                // local to current sampling record s
                Vec3 local_wi = sampling_records[s].uvn.world_to_local(matrix_clusters[slice][c].slice_wi);
                if (local_wi.z() <= 0.0f) continue;

                tree_node_data[num_visible_clusters].p = direction_to_square(local_wi);
                tree_node_data[num_visible_clusters].index = slice * max_slice_clusters + c;    // the radiance data is from (slice, c)
                num_visible_clusters++;
            }
        }                
        
        slices[s].kdtree = new KdTree2(slices[s].tree_node_data, num_visible_clusters, scene);
    }
    stats.toc();
    stats.print_elapsed_milliseconds();

    // TEST:
    // compare slice representative's pixel distribution with the groundtruth distribution at that pixel
    // expect: slice representative (matrix truly evaluated) must resembles the groundtruth!!!
    // the files produced by both of the functions must be exactly the same.
    /*
    int slice = 156;
    int pixel_index = node_data[slices[slice].representative].index;
    int ydebug = pixel_index / 512;
    int xdebug = pixel_index % 512;
    build_map_receiver(ydebug, xdebug);
    build_map(slice);
    */

    Log::info() << "Build Jensen's map..." << endn;
    stats.tic();
    for (int s = 0; s < matrix_clusters.size(); ++s) {
        // build Jensen's map        
        build_map(s);
    }
    stats.toc();
    stats.print_elapsed_milliseconds();
    
    Log::info() << "Estimate sums..." << endn;
    stats.tic();
    // init and estimate sums for each sampling record placed at each slice
    Log::info() << "Density radius: " << density_radius << endn;
    for (int s = 0; s < matrix_clusters.size(); ++s) {
        SamplingRecord &record = sampling_records[s];
        record.start = true;
        record.w = Vec3(0.0f);
        record.f = 0.0f;

        record.sum_f = 0.0f;        // set to zero first
        
        if (radius_type == LocalDensity::NEAREST_NEIGHBOR)
            // too few samples will cause the sum to be inaccurate and patch artifacts can occur.
            record.sum_f = estimate_sum_random(rd, s, num_sum_samples);
        else {
            record.sum_vpl = estimate_sum(s);
           // record.sum_f = estimate_sum_random(rd, s, num_sum_samples);;
           record.sum_f = record.sum_vpl;
        }
        // uvn is initialized when the reduced matrix is established above.

        if (record.sum_f <= 0.0f) {
            record.is_valid = false;
        } else {
            record.is_valid = true;
        }
    }

    stats.toc();
    stats.print_elapsed_milliseconds();

    init_stats.toc();
    init_stats.print_elapsed_seconds();
    Log::info() << "Initialization finished." << endn;
}

Rgb LightSlice::get_slice_color(const Receiver &r) {
    Vec6 p = map_vec6(r.p, r.shading_n, r.m, r.wo, scene_bb_diag);

    int slice;
    //kdtree->find_cluster(p, slice);

    // test find nearest cluster function
    
    vector<int> &slice_indices = metropolis_tmp_slice_indices;
    slice_indices.clear();
    kdtree->find_nearest_cluster(p, metropolis_neighbor_slices, slice_indices);

    slice = slice_indices[0];
    
    
        
    Rgb color;
    if (slice < 0)
        color = DefaultRgb::black;
    else
        color = slices[slice].color;

    return color;    
}

Rgb LightSlice::get_slice_sum(const Receiver &r) {
    Vec6 p = map_vec6(r.p, r.shading_n, r.m, r.wo, scene_bb_diag);

    int slice;
    kdtree->find_cluster(p, slice);

    return sampling_records[slice].sum_f;
}

Rgb LightSlice::gather_clusters(const Receiver &r, IVirtualPointLightEvaluator *vpl) {    
    Vec6 p = map_vec6(r.p, r.shading_n, r.m, r.wo, scene_bb_diag);

    int slice;
    kdtree->find_cluster(p, slice);
    
    if (r.m == NULL) return DefaultRgb::black;
    
    Rgb radiance; 
    Clusters &clusters = matrix_clusters[slice];
    for (int i = 0; i < clusters.size(); ++i) {
        int k = clusters[i].representative;
        
        BrdfPointLight &pl = (*all_vpls)[k];

        bool visible;
        
        Rgb Lo = vpl->radiance(scene, r, pl, visible);
        if (visible) {        
            radiance += clusters[i].weight * Lo;
        }

        //Rgb Lo = vpl->incident_radiance(scene, r, pl);
        //radiance += clusters[i].weight * Lo;
    }

    //if (use_poly_interpolation)
        //fit_polynomial(r, poly_deg, global_incident_cache[slice], coeffs);

    return radiance;
}

/**
 * Sample from slice's representative.
 */
void LightSlice::sample_slice_cluster_cone(int slice, Vec3 &wi, Float &pdf_wi) {

    Clusters &clusters = matrix_clusters[slice];

    Random &rd = *scene->get_random();    
    DiscretePdf &pdf = slices[slice].cluster_pdf;

    int k;
    Float pdf_k;
    pdf.sample(rd, k, pdf_k);

    if (k < 0) {
        wi = Vec3(0.0f);
        pdf_wi = 0.0f;
        return;
    }

    Cluster &cluster = clusters[k];
    Float half_angle = std::min(HALF_PI, std::max(uniform_cone_half_angle, cluster.bc.half_angle()));

    // assume a cone over cluster representative as VSL and perform uniform sampling in the cone
    Vec2 d;
    d.random(rd);
        
    Float cos_theta_max = cos(half_angle);
    
    Float cos_theta = 1.0f - d[0] * (1.0f - cos_theta_max);
    Float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    Float phi = TWO_PI * d[1];
    Float theta = acos(cos_theta);
    
    BrdfPointLight &vpl = (*all_vpls)[cluster.representative];
    Onb uvn;
    uvn.init_from_n(cluster.slice_wi);
    Vec3 cone_dir = uvn.spherical_to_world(theta, phi);

    wi = cone_dir;
    
    pdf_wi = pdf_slice_cluster_cone(slice, wi);
}

/** 
 * Pdf of wi according to sampling from slice representative
 */
Float LightSlice::pdf_slice_cluster_cone(int slice, const Vec3 &wi) {

    DiscretePdf &pdf = slices[slice].cluster_pdf;
    Clusters &clusters = matrix_clusters[slice];  

    Float pdf_wi = 0.0f;
    for (int i = 0; i < clusters.size(); ++i) {        
        if (! clusters[i].is_upper_hemisphere) continue;

        Vec3 dir = clusters[i].slice_wi;
        
        // reuse half angle from cluster cone
        Float half_angle = std::min(HALF_PI, std::max(uniform_cone_half_angle, clusters[i].bc.half_angle()));
        
        Cone cone(dir, half_angle);                
        Float cos_theta_max = cos(half_angle);    
        Float pdf_in_cone = 1.0f / (TWO_PI * (1.0f - cos_theta_max));

        if (cone.contains(wi)) {
            pdf_wi += pdf.probability(i) * pdf_in_cone;
        }
    }
    return pdf_wi;
}

/**
 * In fact, this function implements incoming radiance evaluation. 
 *
 * I(w) is mapped to I(x) where x is in a unit square for evaluation.
 * It is a one-to-one mapping, and I(w) = I(x). No INV_2PI compensation required
 * as there is no integration here.
 */
Float LightSlice::eval_incoming_radiance(int slice, Vec2 &x) {
    //return 1.0f;
    // test: just use Jensen's map
    /*
    int width = luminance_map_resolution;
    int height = luminance_map_resolution;

    int col = (int)(x.x() * width);
    int row = (int)(x.y() * height);
    if (col >= width) col = width - 1;
    if (row >= height) row = height - 1;

    Float pdf_k = slices[slice].luminance_map_pdf.probability(row, col);

    int num_pixels = height * width;

    //Float pixel_area = 1.0f / num_pixels;
    //Float pdf_unit_square = 1.0f / pixel_area * pdf_k;
    Float pdf_unit_square = num_pixels * pdf_k;
    Float pdf_unit_disk = pdf_unit_square * INV_PI;
    Float pdf_wi = pdf_unit_disk * 0.5f;
    return pdf_wi;
    */

    // use 2D-tree to find k-nearest cluster representatives
    KdTree2 *kdtree = slices[slice].kdtree;    
    NodeData2Heap heap(x);

    Float radius = 0.0f;

    if (radius_type == LocalDensity::FIXED_RADIUS) {

        kdtree->find_nearest(x, density_radius, heap);
        
        radius = density_radius;
                
    } else {

        kdtree->find_nearest(x, 4, heap);
        
        if (heap.size() > 0) {
            NodeData2 *furthest = heap.top();
            int s = furthest->index / max_slice_clusters;
            int c = furthest->index % max_slice_clusters;
            Cluster &cluster = matrix_clusters[s][c];

            // local using current slice
            Vec3 local_wi = sampling_records[slice].uvn.world_to_local(cluster.slice_wi);
            radius = (direction_to_square(local_wi) - x).length();
        }
    }
    
    NodeData2Ptrs *nodes = reinterpret_cast<NodeData2Ptrs *>(&heap);

    // take average
    Rgb sum;
    for (int i = 0; i < nodes->size(); ++i) {
        int s, c;
        s = (*nodes)[i]->index / max_slice_clusters;
        c = (*nodes)[i]->index % max_slice_clusters;

        Cluster &cluster = matrix_clusters[s][c];
        if (cluster.is_upper_hemisphere == false) continue;

        Rgb radiance = cluster.weight * R.get_matrix()(s, cluster.representative);
        
        // NOTE: we should not use cluster weight here because this is purely directional radiance interpolation???
        //Rgb radiance = R.get_matrix()(s, cluster.representative);
                
        // distance in the unit square determines weight
        // local using the current slice's basis, not s.
        //Vec3 local_wi = sampling_records[slice].uvn.world_to_local(cluster.slice_wi);
        //Float dist2 = (direction_to_square(local_wi) - x).squared_length();
        //Float weight = exp(dist2 * inv_sigma2);
        
        sum += radiance;
        
    }
    
    if (radius_type == LocalDensity::NEAREST_NEIGHBOR) {
        if (radius > 0.0f)
            sum /= (radius * radius * A_PI);
    }

    // make sure density is not zero so each direction can be selected
    // NOTE: if density_epsilon is used differently, make sure the sum is correctly estimated.
    //Float density = density_epsilon + sum.luminance();    
    //Float density = std::max(density_epsilon_ratio * sampling_records[slice].sum_vpl, sum.luminance());
    //density += density_epsilon;
    Float density = density_epsilon_ratio * sampling_records[slice].sum_vpl + sum.luminance();
    return density;
}

Float LightSlice::estimate_sum(int slice) {
    //return 1.0f;
    // NOTE: this sum estimation must be compatible with the behavior of the eval_incoming_radiance function.
    // This function is only for fixed radius. 

    // we can estimate simply by finding average radiance from the clusters, 
    // and then multiply with the square area, which is 1.    

    Float sum = 0.0f;
    Rgb total_radiance;
              
    int s = slice;
                
    if (! sampling_records[s].is_valid) return 0.0f;

    for (int c = 0; c < matrix_clusters[s].size(); ++c) {

        Cluster &cluster = matrix_clusters[s][c];                

        if (cluster.is_upper_hemisphere) {
            total_radiance += (density_radius * density_radius * A_PI) * cluster.weight * R.get_matrix()(s, cluster.representative);
        }
    }
    

    //sum = density_epsilon * 1.0f /*area */ + total_radiance.luminance();
    sum = total_radiance.luminance();
    return sum * TWO_PI;
}

/**
 * This estimation is better for samples that lies near to the square boundary 
 * (which the disk might be cutoff.)
 */
Float LightSlice::estimate_sum_random(Random &rd, int slice, int num_sum_samples) {

    Float sum = 0.0f;
    for (int i = 0; i < num_sum_samples; ++i) {
        Vec2 s;
        s.random(rd);

        sum += eval_incoming_radiance(slice, s);
    }
    // the TWO_PI compensates for the integral estimated in the unit square instead of its original hemisphere domain
    // we have: dw = 2pi ds or p(w) = p(s) / 2pi. 
    // It is only TWO_PI because we use uniform hemisphere to uniform disk mapping, so p(w) = p(disk) / 2. 
    return sum / num_sum_samples * TWO_PI;
}

Rgb LightSlice::multi_slice_sample(const Receiver &r, Vec3 &wi, Float &pdf_wi) {
    if (!r.m) {
        wi = Vec3(0.0f);
        pdf_wi = 0.0f;
        return DefaultRgb::black;
    }

    // TODO: insert slice representatives into a kd-tree and query instead of doing map and search for slices here. 
    // Caveat: the nearest cache might not have same normals as the point group into a slice.
    Vec6 p = map_vec6(r.p, r.shading_n, r.m, r.wo, scene_bb_diag);
    
    vector<int> &slice_indices = metropolis_tmp_slice_indices;
    slice_indices.clear();
    kdtree->find_nearest_cluster(p, metropolis_neighbor_slices, slice_indices);

    Float min_dist = FLT_MAX;
    int min_slice = -1;
    // TODO: delete slice_pdf here (it's too complicated and too slow if for every pixel we need to look for 
    // a set of neighbors
    vector<Float> &slice_pdf = metropolis_tmp_slice_pdf;
    slice_pdf.resize(slice_indices.size());
    for (int i = 0; i < slice_indices.size(); ++i) {
        if (sampling_records[slice_indices[i]].is_valid && 
            is_view_similar(sampling_records[slice_indices[i]].uvn, r.shading_uvn)) {

            // need to favor the nearest slice presentative (not the parent slice)
            // to minimize distribution difference.
            Vec6 &p = node_data[slices[slice_indices[i]].representative].p; 
            Float dist = (Vec3(p[0], p[1], p[2]) - r.p).length();
            if (dist < min_dist) {
                min_dist = dist;
                min_slice = slice_indices[i];
            }
            slice_pdf[i] = 1.0f / sqrt(distance_epsilon + (Vec3(p[0], p[1], p[2]) - r.p).length());

        } else {

            slice_pdf[i] = 0.0f;
        }
    }

    Random &rd = *scene->get_random();
    /*
    DiscretePdf &pdf = metroplis_tmp_discrete_pdf;
    pdf.set_distribution(slice_pdf);
    
    if (pdf.is_valid() == false) {
    */
    if (min_slice < 0) {
        // there is no suitable slice

        // DEBUG: to see why we have splotches (patchy results)
        /*
        wi = Vec3(0.0f);
        pdf_wi = 0.0f;

        if (sampling_records[slice_indices[0]].is_valid == false)
            return DefaultRgb::red;

        if (is_view_similar(sampling_records[slice_indices[0]].uvn, r.shading_uvn) == false) 
            return DefaultRgb::blue;
        */
        
        // fail due to invalid slice sampling record
        // ask for neighbor slices to help

        vector<int> &slice_indices = metropolis_tmp_slice_indices;
        slice_indices.clear();
        kdtree->find_nearest_cluster(p, metropolis_neighbor_slices_when_failed, slice_indices);

        Float min_dist = FLT_MAX;        
        for (int i = 0; i < slice_indices.size(); ++i) {
            if (sampling_records[slice_indices[i]].is_valid && 
                is_view_similar(sampling_records[slice_indices[i]].uvn, r.shading_uvn)) {

                // pick the nearest valid sampling record to use
                Vec6 &p = node_data[slices[slice_indices[i]].representative].p; 
                Float dist = (Vec3(p[0], p[1], p[2]) - r.p).length();
                if (dist < min_dist) {
                    min_dist = dist;
                    min_slice = slice_indices[i];
                }
            }
        }

        if (min_slice < 0) {
            // still no valid sampling record
            // then no choice, uniform sampling

            Vec2 sampler;
            sampler.random(rd);
            Float phi = TWO_PI * sampler.x();
            Float theta = acos(sampler.y());
            wi = unit_vector(r.shading_uvn.spherical_to_world(theta, phi));
            pdf_wi = INV_2PI;
            return r.m->eval(LocalGeometry(r), r.wo, wi);
        }
    } 

    // good, have a slice to refer to now

    /*
    int k;
    Float pdf_k;
    pdf.sample(rd, k, pdf_k);

    int slice = slice_indices[k];
    */
    // The pdf_k above is accounted in the pdf function.

    int slice = min_slice;
    SamplingRecord &record = sampling_records[slice];
    
    Float pdf_w;
    switch (sampling_type) {
    case LocalSampler::METROPOLIS: 
        metropolis_fetch_cached_sample(slice, wi, pdf_w);            
        break;
            
    case LocalSampler::JENSEN:
        sample_map(slice, wi, pdf_w);
        break;
    }

    if (same_local_wi) {
        wi = r.shading_uvn.local_to_world(record.uvn.world_to_local(wi));    
    }
    
        
    // this direction could have been generated from other slices    
    pdf_wi = multi_slice_pdf(r, wi);

    Rgb brdf = r.m->eval(LocalGeometry(r), r.wo, wi);    
    return brdf;
}


/**
 * Metropolis sampler
 */
void LightSlice::metropolis_fetch_cached_sample(int slice, Vec3 &wi, Float &pdf_wi) {
    
    SamplingRecord &record = sampling_records[slice];

    // assume record is valid and view is similar
    if (! record.is_valid) {

        Log::info() << "This case should not happen. Record and view must have been checked." << endn;
    }

    Random &rd = *scene->get_random();   
    
    const int num_samples = 1024;
    vector<BrdfSample> &samples = slices[slice].samples;

    if (record.start) {
        // run a pass to just warm-up and reduce bias        
        samples.resize(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            metropolis_sample(slice, wi, pdf_wi);
            samples[i] = BrdfSample(wi, pdf_wi, DefaultRgb::black);
                        
            record.start = false;
        }

        // set remains still zero as these samples are just for warming up.
        slices[slice].sample_remains = 0; 
    }

    if (slices[slice].sample_remains == 0) {
        samples.resize(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            metropolis_sample(slice, wi, pdf_wi);
            samples[i] = BrdfSample(wi, pdf_wi, DefaultRgb::black);

            write_to_histogram(slice, wi);
        }

        // random permutation        
        for (int i = 0; i < num_samples; ++i) {
            int k = rd() * num_samples;
            k = std::min(num_samples - 1, k);
            std::swap(samples[i], samples[k]);
        }

        slices[slice].sample_remains = num_samples;
    }

    BrdfSample &sample = samples[slices[slice].sample_remains - 1];
    wi = sample.w;
    pdf_wi = sample.pdf;
    slices[slice].sample_remains--;
}

void LightSlice::metropolis_sample(int slice, Vec3 &wi, Float &pdf_wi) {
    
    Random &rd = *scene->get_random();

    // get the record at the slice for sampling
    SamplingRecord &record = sampling_records[slice];
    
    // assume the record is already valid

    if (record.start) {
        
        // start by choosing a VPL and make a direction
        int k; 
        Float pdf_k;
        slices[slice].cluster_pdf.sample(rd, k, pdf_k);

        if (k < 0) {
            
            Log::info() << "This case should no longer occur." << endn;
            wi = Vec3(0.0f);
            pdf_wi = 0.0f;
            return;

        } else {

            wi = matrix_clusters[slice][k].slice_wi;

        }

        Vec3 local_wi = record.uvn.world_to_local(wi);

        record.start = false;
        record.w = wi;
        record.f = eval_incoming_radiance(slice, direction_to_square(local_wi));
                
        pdf_wi = record.f / record.sum_f;

        return;
    }

    // propose a mutation    
    int mutation_type;
    Float mutation_type_pdf;
    mutation_pdf.sample(rd, mutation_type, mutation_type_pdf);
    
    if (mutation_type == Mutation::UNIFORM) {

        Vec2 sampler;
        sampler.random(rd);
        Float phi = TWO_PI * sampler.x();
        Float theta = acos(sampler.y());
        wi = unit_vector(record.uvn.spherical_to_world(theta, phi));            // all wi calculate relative to slice frame for now
        
    }
    else if (mutation_type == Mutation::VPL) {
        /*
        int k;
        Float pdf_k;
        slices[slice].cluster_pdf.sample(rd, k, pdf_k);

        if (k < 0) {
            // fall back to uniform
            Vec2 sampler;
            sampler.random(rd);
            Float phi = TWO_PI * sampler.x();
            Float theta = acos(sampler.y());
            wi = unit_vector(record.uvn.spherical_to_world(theta, phi));

        } else {
            Cluster &cluster = matrix_clusters[slice][k];

            Vec2 sampler;
            sampler.random(rd);
            Float cos_theta_max = cos(cluster.bc.half_angle());     // use the estimated half angle
            Float cos_theta = 1.0f - sampler.x() * (1.0f - cos_theta_max);
            Float phi = TWO_PI * sampler.y();
            Float theta = acos(cos_theta);

            Onb uvn;
            uvn.init_from_n(cluster.slice_wi);
        
            wi = uvn.spherical_to_world(theta, phi);
        }
        */

        sample_slice_cluster_cone(slice, wi, pdf_wi);

    }
    else {     // perturb

        Vec2 sampler;
        sampler.random(rd);
        Float cos_theta = 1.0f - sampler.x() * (1.0f - perturb_cos_theta_max);
        Float phi = TWO_PI * sampler.y();
        Float theta = acos(cos_theta);

        Onb uvn;
        uvn.init_from_n(record.w);
        
        wi = uvn.spherical_to_world(theta, phi);

    }
    
    Float f = 0.0f;
        
    Vec3 local_wi = record.uvn.world_to_local(wi);
    if (local_wi.z() > 0.0f) {    // prevent jumping to lower half of the hemisphere
    
        Vec2 x = direction_to_square(local_wi);
        f = eval_incoming_radiance(slice, x);
                
        Float transition_to_new = 1.0f;
        Float transition_to_old = 1.0f;
        
        transition_to_new = mutation_pdf.probability(Mutation::UNIFORM) * INV_2PI + 
                            mutation_pdf.probability(Mutation::CONE_PERTURB) * perturb_cone_pdf +
                            (mutation_pdf.probability(Mutation::VPL) > 0.0f ? 
                                mutation_pdf.probability(Mutation::VPL) * pdf_slice_cluster_cone(slice, wi) : 0.0f);

        transition_to_old = mutation_pdf.probability(Mutation::UNIFORM) * INV_2PI + 
                            mutation_pdf.probability(Mutation::CONE_PERTURB) * perturb_cone_pdf +
                            (mutation_pdf.probability(Mutation::VPL) > 0.0f ? 
                                mutation_pdf.probability(Mutation::VPL) * pdf_slice_cluster_cone(slice, record.w) : 0.0f);
                                
        Float acceptance = std::min(1.0f, (f * transition_to_old) / (record.f * transition_to_new));
        
        Float t = rd();
        if (t < acceptance) {
        
            // save state
            record.w = wi;
            record.f = f;

        } else {    
        
            // revert to state saved in record
            wi = record.w;
            f = record.f;
        }

    } else {

        // revert to state saved in record
        wi = record.w;
        f = record.f;
    }
    
    pdf_wi = f / record.sum_f;
}

Float LightSlice::metropolis_get_pdf(int slice, const Vec3 &wi) {
    
    SamplingRecord &record = sampling_records[slice];
    
    // assume record is valid and view is similar

    // assume same world direction: change to local of slice    
    Vec3 warp_local_wi = record.uvn.world_to_local(wi);
    if (warp_local_wi.z() <= 0.0f) return 0.0f;

    Vec2 x = direction_to_square(warp_local_wi);
    Float f = eval_incoming_radiance(slice, x);

    Float pdf_wi = f / record.sum_f;
    return pdf_wi;
}

Float LightSlice::multi_slice_pdf(const Receiver &r, const Vec3 &wi) {
    if (!r.m) {
        return 0.0f;
    }
            
    // only opaque surface
    if (dot(r.shading_n, wi) <= 0.0f) return 0.0f;
        
    Vec6 p = map_vec6(r.p, r.shading_n, r.m, r.wo, scene_bb_diag);

    vector<int> &slice_indices = metropolis_tmp_slice_indices;
    slice_indices.clear();
    kdtree->find_nearest_cluster(p, metropolis_neighbor_slices, slice_indices);

    vector<Float> &slice_pdf = metropolis_tmp_slice_pdf;
    slice_pdf.resize(slice_indices.size());
    
    Float min_dist = FLT_MAX;
    int min_slice = -1;

    bool is_zero = true;;
    for (int i = 0; i < slice_indices.size(); ++i) {
        if (sampling_records[slice_indices[i]].is_valid && 
            is_view_similar(sampling_records[slice_indices[i]].uvn, r.shading_uvn)) {

            Vec6 &p = node_data[slices[slice_indices[i]].representative].p;    
            
            Float dist = (Vec3(p[0], p[1], p[2]) - r.p).length();
            if (dist < min_dist) {
                min_dist = dist;
                min_slice = slice_indices[i];
            }
            slice_pdf[i] = 1.0f / sqrt(distance_epsilon + (Vec3(p[0], p[1], p[2]) - r.p).length());

            is_zero = false;
        } else {

            slice_pdf[i] = 0.0f;
        }
    }

    if (min_slice < 0) {
        vector<int> &slice_indices = metropolis_tmp_slice_indices;
        slice_indices.clear();
        kdtree->find_nearest_cluster(p, metropolis_neighbor_slices_when_failed, slice_indices);

        Float min_dist = FLT_MAX;        
        for (int i = 0; i < slice_indices.size(); ++i) {
            if (sampling_records[slice_indices[i]].is_valid && 
                is_view_similar(sampling_records[slice_indices[i]].uvn, r.shading_uvn)) {

                // pick the nearest valid sampling record to use
                Vec6 &p = node_data[slices[slice_indices[i]].representative].p; 
                Float dist = (Vec3(p[0], p[1], p[2]) - r.p).length();
                if (dist < min_dist) {
                    min_dist = dist;
                    min_slice = slice_indices[i];
                }
            }
        }

        if (min_slice < 0) return INV_2PI;

    }


    //if (is_zero) return INV_2PI;
    //if (is_zero) return 0.0f;
    //if (is_zero) return pdf_cluster_cone(r, wi);
    
    // TEST: use min
    SamplingRecord &record = sampling_records[min_slice];
    Vec3 dir;
    if (same_local_wi) {
        dir = record.uvn.local_to_world(r.shading_uvn.world_to_local(wi));    // warped
    }
    else {
        dir = wi;
    }
    switch (sampling_type) {
    case LocalSampler::METROPOLIS:
        return metropolis_get_pdf(min_slice, dir);
        
    case LocalSampler::JENSEN:
        return pdf_map(min_slice, dir);
    }
    // end of function here

    DiscretePdf &pdf = metroplis_tmp_discrete_pdf;
    pdf.set_distribution(slice_pdf);
    
    Float pdf_wi = 0.0f;
    for (int i = 0; i < slice_indices.size(); ++i) {
        
        int slice = slice_indices[i];
        SamplingRecord &record = sampling_records[slice];
        
        // wi could have been generated by this slice

        if (sampling_records[slice_indices[i]].is_valid && 
            is_view_similar(sampling_records[slice_indices[i]].uvn, r.shading_uvn)) {

            // no matter whether the assumption is to have same world or local wi,
            // we still must ensure that distributions at triangle boundaries to be the same, i.e., 
            // boundary cache points must have same local basis. Otherwise their distributions are still rotated,
            // and artifacts will occur.
            Vec3 dir;
            if (same_local_wi) {
                dir = record.uvn.local_to_world(r.shading_uvn.world_to_local(wi));    // warped
            } else {
                dir = wi;
            }
            
            switch (sampling_type) {
            case LocalSampler::METROPOLIS: 
                pdf_wi += pdf.probability(i) * metropolis_get_pdf(slice, dir);
                break;
            case LocalSampler::JENSEN:
                pdf_wi += pdf.probability(i) * pdf_map(slice, dir);
                break;
            }
        }
    }
    
    return pdf_wi;
}


void LightSlice::build_map(int slice) {    
    ImageFloat &luminance_map = slices[slice].luminance_map;        
    luminance_map.initialize(0.0f, luminance_map_resolution, luminance_map_resolution, 1);
    
    int width = luminance_map_resolution;
    int height = luminance_map_resolution;
    
    Float sum = 0.0f;
    int s = slice;
    Clusters &clusters = matrix_clusters[s];
                
    if (sampling_records[s].is_valid) {
        
        for (int c = 0; c < clusters.size(); ++c) {
            if (! clusters[c].is_upper_hemisphere) continue;
            //equiv to (visibility > 0.0f && dot(recv.shading_n, wi) > 0.0f && incident_power != DefaultRgb::black) 

            Vec3 wi = clusters[c].slice_wi;

            // local to current slice, not s
            Vec3 local_wi = sampling_records[slice].uvn.world_to_local(wi);
            if (local_wi.z() <= 0.0f) continue;
                        
            // but retrieve incident radiance to s
            Rgb radiance = clusters[c].weight * R.get_matrix()(s, clusters[c].representative);

            Vec2 sq = direction_to_square(local_wi);
            int x = (int)(sq.x() * width);
            int y = (int)(sq.y() * height);
            if (x >= width) x = width - 1;
            if (y >= height) y = height - 1;
            if (x < 0) {
                Log::info() << "Unexpected case." << endn;
                x = 0;
            }
            if (y < 0) {
                Log::info() << "Unexpected case." << endn;
                y = 0;
            }
            luminance_map(y, x) += radiance.luminance();
            sum += radiance.luminance();
        
            // 
            // ARTIFACT test: 
            /*
            Vec3 wi = clusters[c].slice_wi;
            Vec3 local_wi = sampling_records[slice].uvn.world_to_local(wi);
            if (local_wi.z() <= 0.0f) continue;

            // but retrieve incident radiance to s
            Rgb radiance = clusters[c].weight * R.get_matrix()(s, clusters[c].representative);

            Vec2 sq = direction_to_square(local_wi);
            Float pdf = eval_incoming_radiance(s, sq);
            
            int x = (int)(sq.x() * width);
            int y = (int)(sq.y() * height);
            if (x >= width) x = width - 1;
            if (y >= height) y = height - 1;
            if (x < 0) {
                Log::info() << "Unexpected case." << endn;
                x = 0;
            }
            if (y < 0) {
                Log::info() << "Unexpected case." << endn;
                y = 0;
            }
            luminance_map(y, x) += pdf;
            sum += pdf;
            */
        }
    

        // avoid bias by setting some epsilon to zero cells
        for (int j = 0; j < width; ++j) {
            for (int i = 0; i < height; ++i) {
                // this can still fail if the sum is zero.
                luminance_map(i, j) += luminance_map_epsilon * sum;
            }
        }

    }  // end if

    ostringstream name2;
    name2 << "data/luminance_map_" << slice << ".exr";
    //luminance_map.save(name2.str().c_str());

    slices[slice].luminance_map_pdf.set_distribution(luminance_map);

    // initialize histogram
    ImageFloat &histogram = slices[slice].histogram;
    histogram.initialize(0.0f, luminance_map_resolution, luminance_map_resolution, 1);
    frame_index = 0;
}

void LightSlice::save_map(int slice) {
    ImageFloat color_map;
    
    color_map.initialize(0.0f, luminance_map_resolution, luminance_map_resolution, 3);
        
    int width = color_map.get_size().width;
    int height = color_map.get_size().height;
    
    int s = slice;
    Clusters &clusters = matrix_clusters[s];
                
    if (! sampling_records[s].is_valid) return;
        
    for (int c = 0; c < clusters.size(); ++c) {
        if (! clusters[c].is_upper_hemisphere) continue;
            
        Vec3 wi = clusters[c].slice_wi;

        // local to current slice, not s
        Vec3 local_wi = sampling_records[slice].uvn.world_to_local(wi);
        if (local_wi.z() <= 0.0f) continue;
                        
        // but retrieve incident radiance to s
        Rgb radiance = clusters[c].weight * R.get_matrix()(s, clusters[c].representative);

        Vec2 sq = direction_to_square(local_wi);
        int x = (int)(sq.x() * width);
        int y = (int)(sq.y() * height);
        if (x >= width) x = width - 1;
        if (y >= height) y = height - 1;
                        
        color_map(y, x, 0) += radiance.red();
        color_map(y, x, 1) += radiance.green();
        color_map(y, x, 2) += radiance.blue();
    }
    

    ostringstream name;
    name << "data/color_map_" << slice << ".exr";
    color_map.save(name.str().c_str());

    ImageFloat interpolate_luminance;
    interpolate_luminance.initialize(0.0f, height, width, 1);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            Float x = (j + 0.5f) / width;
            Float y = (i + 0.5f) / height;
            Vec2 s(x, y);
            interpolate_luminance(i, j) = eval_incoming_radiance(slice, s);
        }
    }
    
    ostringstream name3;
    name3 << "data/interpolate_map_" << slice << ".exr";
    interpolate_luminance.save(name3.str().c_str());
}

/**
 * Build the density maps using VPLs to a receiver.
 *
 * For comparison.
 */
void LightSlice::build_map_receiver(int row, int col) {
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float tick = scene->get_tick();
    Float max_tmax = scene->get_max_tmax();
    Camera *camera = scene->get_camera();
    Ray ray = camera->shoot_ray(row, col);
    HitRecord rec;
    if (!agg->hit(ray, tmin, max_tmax, tick, rec) || rec.material == NULL) {
        return;
    }

    Receiver recv(ray, rec);
    Vec6 p = map_vec6(recv.p, recv.shading_n, recv.m, recv.wo, scene_bb_diag);
    
    int slice;
    kdtree->find_cluster(p, slice);

    Clusters &clusters = matrix_clusters[slice];

    // check for view difference
    Onb &uvn = sampling_records[slice].uvn;
    bool is_similar = is_view_similar(uvn, recv.shading_uvn);
    
    ImageFloat receiver_map;
    receiver_map.initialize(0.0f, luminance_map_resolution, luminance_map_resolution, 3);

    int width = luminance_map_resolution;
    int height = luminance_map_resolution;

    for (int c = 0; c < clusters.size(); ++c) {
        BrdfPointLight &vpl = (*all_vpls)[clusters[c].representative];
        if (agg->hit(recv.p, vpl.org(), tmin, tick)) continue;
        /*
        // only for verifying slice's visibility
        Float visible = matV(slice, clusters[c].representative);
        if (visible == 0.0f) {
            Log::info() << "Assert visibility failed. Visibility should be 1." << endn;
        }*/

        Vec3 wi = unit_vector(vpl.org() - recv.p);
        //if (recv.m->pdf(LocalGeometry(recv), recv.wo, wi) <= 0.0f) continue;
        // we don't care about r.wo, so do not use pdf which might check if it is under surface
        if (dot(recv.shading_n, wi) <= ZERO_EPSILON) continue;

        Vec3 local_wi = recv.shading_uvn.world_to_local(wi);
        Vec2 s = direction_to_square(local_wi);
        int x = (int)(s.x() * width);
        int y = (int)(s.y() * height);
        if (x >= width) x = width - 1;
        if (y >= height) y = height - 1;
        
        // NOTE: should we use cluster.weight here?        
        Rgb Lin = evaluator->incident_radiance(scene, recv, vpl);
        Rgb radiance = clusters[c].weight * Lin;
        
        receiver_map(y, x, 0) += radiance.red();
        receiver_map(y, x, 1) += radiance.green();
        receiver_map(y, x, 2) += radiance.blue();
    }    

    
    ostringstream name;
    name << "data/receiver" << "_row" << row << "_col" << col << "_slice" << slice;
    if (is_similar == false)
        name << "_viewdiff";
    name << ".exr";
    receiver_map.save(name.str().c_str());
    
}

void LightSlice::write_to_histogram(int slice, Vec3 &wi) {
    ImageFloat &histogram = slices[slice].histogram;
    
    Vec3 local_wi = sampling_records[slice].uvn.world_to_local(wi);
    Vec2 s = direction_to_square(local_wi);
    int width = luminance_map_resolution;
    int height = luminance_map_resolution;

    int col = (int)(s.x() * width);
    int row = (int)(s.y() * height);
    if (col >= width) col = width - 1;
    if (row >= height) row = height - 1;

    histogram(row, col) += 1.0f;
}

void LightSlice::save_histogram(int slice) {
    ostringstream name;
    name << "data/histogram_" << slice << ".exr";
    slices[slice].histogram.save(name.str().c_str());
}

void LightSlice::on_frame_end() {
    
    
}

void LightSlice::sample_map(int slice, Vec3 &wi, Float &pdf_wi) {
    Random &rd = *scene->get_random();

    if (slices[slice].luminance_map_pdf.is_valid() == false) {
        Vec2 sampler;
        sampler.random(rd);
        Float phi = TWO_PI * sampler.x();
        Float theta = acos(sampler.y());
        wi = unit_vector(sampling_records[slice].uvn.spherical_to_world(theta, phi));
        pdf_wi = INV_2PI;
        return;
    }

    int row, col;
    Float pdf_k;
    slices[slice].luminance_map_pdf.sample(rd, row, col, pdf_k);

    if (row < 0 || col < 0) {
        Log::info() << "This case should no longer occur!" << endn;
        wi = Vec3(0.0f);
        pdf_wi = 0.0f;
        return;
    }

    Vec2 jitter;
    jitter.random(rd);

    Vec2 s((col + jitter.x()) / luminance_map_resolution, (row + jitter.y()) / luminance_map_resolution);
    Vec3 local_wi = square_to_direction(s);
    wi = sampling_records[slice].uvn.local_to_world(local_wi);

    int num_pixels = luminance_map_resolution * luminance_map_resolution;
    Float pdf_unit_square = num_pixels * pdf_k;
    Float pdf_unit_disk = pdf_unit_square * INV_PI;
    pdf_wi = pdf_unit_disk * 0.5f;   // total is 1 / 2PI, not 1 / 2PI^2.
}

Float LightSlice::pdf_map(int slice, const Vec3 &wi) {
    Vec3 local_wi = sampling_records[slice].uvn.world_to_local(wi);
    if (local_wi.z() <= 0.0f) return 0.0f;

    if (slices[slice].luminance_map_pdf.is_valid() == false) {
        return INV_2PI;
    }

    Vec2 s = direction_to_square(local_wi);
    int width = luminance_map_resolution;
    int height = luminance_map_resolution;

    int col = (int)(s.x() * width);
    int row = (int)(s.y() * height);
    if (col >= width) col = width - 1;
    if (row >= height) row = height - 1;

    Float pdf_k = slices[slice].luminance_map_pdf.probability(row, col);

    int num_pixels = height * width;

    //Float pixel_area = 1.0f / num_pixels;
    //Float pdf_unit_square = 1.0f / pixel_area * pdf_k;
    Float pdf_unit_square = num_pixels * pdf_k;
    Float pdf_unit_disk = pdf_unit_square * INV_PI;
    Float pdf_wi = pdf_unit_disk * 0.5f;
    return pdf_wi;
}

Rgb LightSlice::sample_uniform(const Receiver &r, Vec3 &wi, Float &pdf_wi) {
    if (! r.m) {
        wi = Vec3(0.0f);
        pdf_wi = 0.0f;
        return DefaultRgb::black;
    }
    
    Random &rd = *scene->get_random();

    Vec2 sampler;
    sampler.random(rd);
    Float phi = TWO_PI * sampler.x();
    Float theta = acos(sampler.y());
    wi = unit_vector(r.shading_uvn.spherical_to_world(theta, phi));
    pdf_wi = INV_2PI;
    return r.m->eval(LocalGeometry(r), r.wo, wi);        
}

Float LightSlice::pdf_uniform(const Receiver &r, const Vec3 &wi) {
    if (! r.m) return 0.0f;
    if (dot(r.shading_n, wi) <= 0.0f) return 0.0f;

    return INV_2PI;
}

Rgb LightSlice::sample(Random &rd, const Receiver &r, Vec3 &wi, Float &pdf_wi) {
    return multi_slice_sample(r, wi, pdf_wi);
    //return sample_cluster_cone(r, wi, pdf_wi);
    //return sample_uniform(r, wi, pdf_wi);
}

Float LightSlice::pdf(const Receiver &r, const Vec3 &wi) {
    return multi_slice_pdf(r, wi);
    //return pdf_cluster_cone(r, wi);
    //return pdf_uniform(r, wi);    
}

bool LightSlice::is_valid() {
    return true;
}

bool LightSlice::is_singular() {
    return false;
}

void LightSlice::slice_info(int y, int x) {    
    Aggregate *agg = scene->get_aggregate();
    Float tmin = scene->get_tmin();
    Float tick = scene->get_tick();
    Float max_tmax = scene->get_max_tmax();
    Camera *camera = scene->get_camera();
    Ray ray = camera->shoot_ray(y, x);
    HitRecord rec;
    if (!agg->hit(ray, tmin, max_tmax, tick, rec) || rec.material == NULL) {
        return;
    }

    Receiver recv(ray, rec);
    Vec6 p = map_vec6(recv.p, recv.shading_n, recv.m, recv.wo, scene_bb_diag);

    int slice;
    kdtree->find_cluster(p, slice);

    Log::info() << "Slice: " << slice << endn;
    
    save_map(slice);
    save_histogram(slice);
}

void LightSlice::on_click(MouseButton button, MouseState state, int x, int y) {
    slice_info(y, x);

    build_map_receiver(y, x);
}

} // end namespace