#ifndef _MUTABLE_PATH_H_
#define _MUTABLE_PATH_H_

#include "common.h"
#include "form_factor.h"
#include "local_geometry.h"
#include "material.h"
#include "area_light.h"
#include "scene.h"

namespace Renzoku {

struct SamplingRecord {
};

/**
 * This version of path vertex is different from PathNode used in bidirectional path tracing. 
 *
 * Each vertex caches wi, wo, and radiometric terms so that probability of the path and MIS weight
 * can be easily calculated.
 */
struct PathVertex {    
    LocalGeometry dg;

    Material *material;
    Light    *light;
    
    enum Type {                     // the way this vertex is sampled (so given a vertex we know which value out of the two to use)
        EYE,
        LIGHT
    } type;

    
    // cache values at each vertex
    Rgb throughput;                 // when path is complete, throughput from eye or light is the same

    Rgb importance;                 // throughput / pdf of path sampled from eye                        
    Rgb flux;                       // throughput / pdf of path sampled from light

    Vec3 wi;                        // direction towards next vertex in the subpath
    Vec3 wo;                        // direction towards previous vertex in the subpath
                        
    Float pdf_by_next;              // pdf (area measure) of how this vertex is sampled at first or from it's next vertex.
    Float pdf_by_prev;              // pdf (area measure) of how this vertex is sampled at first or from it's previous vertex.
    Float pdf_dA;                   // the actual area measure

    SamplingRecord sr;              // record sampling decision made to at the vertex

    PathVertex() : material(NULL), light(NULL), type(EYE), pdf_by_prev(0.0f), pdf_by_next(0.0f), pdf_dA(0.0f) {}

    /** 
     * Eye node
     */
    PathVertex(const Vec3 &p, const Vec3 &n, const Rgb &We) : 
        
        material(NULL), light(NULL),
        type(EYE),
        pdf_by_prev(1.0f), pdf_by_next(0.0f), pdf_dA(1.0f)
    {
        
        importance = We;
        flux = DefaultRgb::black;
        throughput = DefaultRgb::black;
        dg.p = p;
        dg.n = n;
        
    }

    /**
     * General surface node.
     */
    PathVertex(Type type, const LocalGeometry &dg, Material *m) :
        dg(dg), 
        material(m), light(NULL),
        type(type), 
        pdf_by_prev(0.0f), pdf_by_next(0.0f), pdf_dA(0.0f)
    {
    }

    /**
     * Light node (starting a light subpath or ending an eye subpath).
     *
     * pdf_light: probablity of the light being selected
     */
    PathVertex(Type type, const Vec3 &p, const Vec3 &n, Light *light, Float pdf_light = 1.0f) :
      
        material(NULL), light(light), 
        type(type),
        pdf_by_prev(pdf_light), pdf_by_next(0.0f), pdf_dA(pdf_light)
    {

        dg.p = p;
        dg.n = n;

    }

    inline bool is_on_light() const {
        return (light != NULL && material == NULL);
    }
    
    inline bool is_on_sensor() const {
        return (light == NULL && material == NULL);
    }

    inline bool is_on_surface() const {
        return (light == NULL && material != NULL);
    }

    bool has_next_vertex() const {
        return wi != Vec3(0.0f);
    }

    bool has_prev_vertex() const {
        return wo != Vec3(0.0f);
    }

};

enum PathVertexConnectionType {
    CONNECTION,
    MERGING
};

const int MAX_MUTABLE_PATH_NODE = 64;

static Float pdf_w_to_area(Float pdf_w, const Vec3 &p, const Vec3 &q, const Vec3 &qn) {
    Vec3 d = p - q;
    return pdf_w * fabs(dot(unit_vector(d), qn)) / (d.squared_length());
}

// compute area probabilty of current vertex sampled from previous vertex
static Float pdf_area(const PathVertex &prev, const PathVertex &current) {
    Float pdf_w;
    if (prev.is_on_light()) {
        if (prev.light->get_light_type() == Light::AREA_LIGHT) {
            AreaLight *area_light = static_cast<AreaLight *>(prev.light);
            pdf_w = area_light->pdf_wo(prev.dg.p, prev.wi);
        }
        else {
            pdf_w = 0.0f;     // not yet supported
        }
    }
    else if (prev.is_on_sensor()) {
        
        //pdf_w = 1.0f;
        throw std::exception("this case is not handled in this function");

    }
    else {
        pdf_w = prev.material->get_bsdf()->pdf(prev.dg.uvn, prev.wo, prev.wi);
    }
    Vec3 d = prev.dg.p - current.dg.p;
    return pdf_w * fabs(dot(unit_vector(d), current.dg.n)) / (d.squared_length());
}

struct MutablePath {    
    PathVertex nodes[MAX_MUTABLE_PATH_NODE];    /// if full path, vertex 0 is on light, and vertex k is on eye. 
    int num_nodes;
    
    int connection_index;                       /// if it is bidirectional, this is where the connection happens
                                                /// vertices involved are i, i + 1.
                                                /// -1: unidirectional light path
                                                /// num_nodes: unidirectional eye path

    
    PathVertexConnectionType connection_type;   /// vertex connection or merging
                                        
    Scene *scene;
    Vec2 pixel; 
    Float pdf_pixel;                    /// area measure of the pixel

    bool is_valid;                      /// when path is marked invalid, it is discarded.
    bool is_complete;                   /// indicates a full light path (with eye and light vertex found)


    /**
     * Create a path where its node buffer is provided from node pool.
     */    
    MutablePath(Scene *scene) : scene(scene) {
		clear();
    }

    ~MutablePath() {
        
    }

    MutablePath(const MutablePath &p) {
        this->operator=(p);
    }

    const MutablePath *operator=(const MutablePath &p) {
        this->scene = p.scene;
        memcpy(this->nodes, p.nodes, sizeof(PathVertex) * MAX_MUTABLE_PATH_NODE);
        this->num_nodes = p.num_nodes;
        this->is_valid = p.is_valid;
        this->is_complete = p.is_complete;
        this->pixel = p.pixel;
        this->pdf_pixel = p.pdf_pixel;
        this->connection_index = p.connection_index;
        return this;
    }

	void clear() {
        num_nodes = 0;
		is_valid = true;
		is_complete = false;
        pixel = Vec2(0.0f);
        pdf_pixel = 0.0f;
        connection_index = -1;
	}

    void set_pixel(const Vec2 &pixel, Float pdf_pixel) {
        this->pixel = pixel;
        this->pdf_pixel = pdf_pixel;
    }
    
    bool _extend(const PathVertex &new_node, bool connect) {

        if (num_nodes >= MAX_MUTABLE_PATH_NODE) return false;

        if (is_complete) return false;

        nodes[num_nodes++] = new_node;

        if (num_nodes == 1) {
            PathVertex &current = nodes[num_nodes - 1];
            if (current.is_on_light()) {
                if (current.light->get_light_type() == Light::AREA_LIGHT) {
                    AreaLight *area_light = static_cast<AreaLight *>(current.light);
                    current.pdf_by_prev *= area_light->pdf(current.dg.p);
                } else {
                    current.pdf_by_prev *= 0.0f;     // not yet supported
                }
            } else if (current.is_on_sensor()) {
                
                // already set in constructor
                // current.pdf_by_prev = 1.0f;

            } else {
                num_nodes--;
                return false;
            }

            if (!connect) {
                current.pdf_dA = current.pdf_by_prev;
            }
            else {
                current.pdf_dA = new_node.pdf_dA;
            }
            return true;
        }

        if (num_nodes == 2) {
            PathVertex &current = nodes[num_nodes - 1];
            PathVertex &prev = nodes[num_nodes - 2];

            current.wo = unit_vector(prev.dg.p - current.dg.p);
            prev.wi = -current.wo;

            // for eye path current should be 1?
            if (prev.is_on_light()) {
                current.pdf_by_prev = pdf_area(prev, current);
                if (current.is_on_sensor()) {    // if current is a surface point, not enough information to sample prev yet (BRDF at current not evaluable).
                    prev.pdf_by_next = pdf_area(current, prev);
                }
            }
            else if (prev.is_on_sensor()) {
                /*
                Camera *camera = scene->get_camera();
                Float focal = camera->get_focal_length();
                Vec3 sensor_point = camera->get_sensor_point(pixel);
                Float cos_alpha = focal / sensor_point.length();
                Float pdf_w = pdf_pixel * focal * focal / (cos_alpha * cos_alpha * cos_alpha);
                current.pdf_by_prev = pdf_w_to_area(pdf_w, prev.dg.p, current.dg.p, current.dg.n);
                */

                Camera *camera = scene->get_camera();
                Float focal = camera->get_focal_length();
                Vec3 sensor_point = camera->get_sensor_point(pixel);
                Float cos_alpha = focal / sensor_point.length();
                Float G = FormFactor::form_factor_abs(prev.dg.p, prev.dg.n, current.dg.p, current.dg.n);
                current.pdf_by_prev = pdf_pixel * focal * focal / pow(cos_alpha, 4) * G;

                if (current.is_on_light()) {    // if current is a surface point, not enough information to sample prev yet (BRDF at current not evaluable).
                    prev.pdf_by_next = pdf_area(current, prev);
                }
            }
            else {
                throw std::exception("not yet handled");
            }
            if (!connect) {
                current.pdf_dA = current.pdf_by_prev;
            }
            else {
                current.pdf_dA = new_node.pdf_dA;
            }

            // manage flux or importance for path of length 1
            if (prev.is_on_light()) {
                Rgb Le;
                if (prev.light->get_light_type() == Light::AREA_LIGHT) {
                    AreaLight *area_light = static_cast<AreaLight *>(prev.light);
                    area_light->query(prev.dg.p, prev.dg.n, prev.wi /*out dir*/, Le);
                }
                else {
                    // not yet supported
                }
                Float G = FormFactor::form_factor_abs(prev.dg.p, prev.dg.n, current.dg.p, current.dg.n);
                Float pdf_area = prev.pdf_dA * current.pdf_dA;
                if (pdf_area > 0.0f) {
                    current.flux = Le * G / pdf_area;
                    current.throughput = Le * G;
                } else {
                    current.flux = DefaultRgb::black;
                    current.throughput = DefaultRgb::black;
                }

                if (current.is_on_sensor())
                    is_complete = true;
            }
            else if (prev.is_on_sensor()) {
                Camera *camera = scene->get_camera();
                Float focal = camera->get_focal_length();
                Vec3 sensor_point = camera->get_sensor_point(pixel);
                Float cos_alpha = focal / sensor_point.length();
                Float image_height = camera->get_image_size().height;
                
                //Float tan_half_fov = tan(camera->get_vertical_fov() * 0.5f * DEGREE_TO_RADIAN);
                //Float pixel_area = 4 * focal * focal * tan_half_fov * tan_half_fov / (image_height * image_height);
                // faster
                //Float film_height = camera->get_sensor()->get_film_size().height;
                //Float pixel_area = (film_height * film_height) / (image_height * image_height);
                
                Float h_weight = 1.0f;    // assume box filter in image target
                Rgb We = focal * focal / pow(cos_alpha, 4) * h_weight;

                Float G = FormFactor::form_factor_abs(prev.dg.p, prev.dg.n, current.dg.p, current.dg.n);
                Float pdf_area = prev.pdf_dA * current.pdf_dA;
                if (pdf_area > 0.0f) {
                    current.importance = We * G / pdf_area;
                    current.throughput = We * G;
                } else {
                    current.importance = DefaultRgb::black;
                    current.throughput = DefaultRgb::black;
                }
                if (current.is_on_light())
                    is_complete = true;
            }
            else {
                
                throw std::exception("should not happen");

            }
            return true;
        }

        // stand at the vertex before the new node
        PathVertex &current = nodes[num_nodes - 2];     
        PathVertex &prev = nodes[num_nodes - 3];
        PathVertex &next = nodes[num_nodes - 1];

        // maintain directions
        //current.wo = unit_vector(prev.dg.p - current.dg.p);
        //prev.wi = -current.wo;
        current.wi = unit_vector(next.dg.p - current.dg.p);
        next.wo = -current.wi;

        // compute its probability pdf_by_prev and pdf_by_next
        /*
        if (prev.is_on_light() || prev.is_on_sensor()) {
            // already computed
        } else {
            current.pdf_by_prev = pdf_area(prev, current);    // already computed????
        }
        current.pdf_by_next = pdf_area(next, current);      // only correct if next is on light or sensor
        // if material, will be corrected when another new node is added
        */

        prev.pdf_by_next = pdf_area(current, prev);
        next.pdf_by_prev = pdf_area(current, next);

        if (next.is_on_light() || next.is_on_sensor()) {
            current.pdf_by_next = pdf_area(next, current);
        }
        


        if (!connect) {
            next.pdf_dA = next.pdf_by_prev;
        }
        else {
            next.pdf_dA = new_node.pdf_dA;
        }

        // compute importance, flux in the last segment
        Float G = 0.0f;
        Rgb delta = DefaultRgb::black;
        Rgb delta_throughput = DefaultRgb::black;

        if (next.pdf_dA > 0.0f) {
            G = FormFactor::form_factor_abs(current.dg.p, current.dg.n, next.dg.p, next.dg.n);
            delta = current.material->get_bsdf()->eval(current.dg.uvn, current.wo, current.wi);
            delta_throughput = delta;

            if (next.is_on_light()) {
                Rgb Le;
                if (next.light->get_light_type() == Light::AREA_LIGHT) {
                    AreaLight *area_light = static_cast<AreaLight *>(next.light);
                    area_light->query(next.dg.p, next.dg.n, next.wo /* hit light */, Le);
                }
                else {
                    // not yet supported
                }
                delta *= Le * G / next.pdf_dA;
                delta_throughput *= Le * G;
            }
            else if (next.is_on_sensor()) {
                
                // TODO: this is not white?
                Rgb We = DefaultRgb::white;

                delta *= We * G / next.pdf_dA;
                delta_throughput *= We * G;
            }
            else {
                delta *= G / next.pdf_dA;
                delta_throughput *= G;
            }   
        
        }

        PathVertex &root = nodes[0];
        if (root.is_on_light()) {
            next.flux = current.flux * delta;
        }
        else {
            next.importance = current.importance * delta;
        }
        next.throughput = current.throughput * delta_throughput;

        if (num_nodes >= 3 &&
            ((root.is_on_light() && next.is_on_sensor()) ||
             (root.is_on_sensor() && next.is_on_light())
            )) {

            is_complete = true;
        }

        return true;
    }


    /**
    * Extend a unidirectional subpath with a new node
    * sampled from the current last node.
    */
    bool extend(const PathVertex &new_node) {
        return _extend(new_node, false);
    }

    /**
     * Connect the current last node of the path to the new node. 
     * The new node is not needed to be sampled from the current last node. 
     */
    bool connect(const PathVertex &new_node) {
        return _extend(new_node, true);
    }

    bool is_light_subpath() {
        return (num_nodes > 0 && nodes[0].is_on_light());
    }

    bool is_eye_subpath() {
        return (num_nodes > 0 && nodes[0].is_on_sensor());
    }

    const PathVertex *get_subpath_last_vertex() const {
        if (num_nodes > 0)
            return &nodes[num_nodes - 1];
        return NULL;
    }

    bool is_unidirectional() {
        return (connection_index == -1);
    }

    /*
    Rgb brdf(int node_index, const Vec3 &wi) {
        if (num_nodes <= 1) return DefaultRgb::black;
        if (node_index < 0 || node_index >= num_nodes) return DefaultRgb::black;

        if (nodes[node_index].material == NULL) return DefaultRgb::black;

        return nodes[node_index].material->eval(nodes[node_index].dg, nodes[node_index].wo, wi);
    }*/

    Rgb brdf(int node_index) {
        // except first and last node which has no material
        if (num_nodes <= 1 || node_index >= num_nodes - 1) return DefaultRgb::black;
        
        return nodes[node_index].material->eval(nodes[node_index].dg, 
                                                nodes[node_index].wo, nodes[node_index].wi);
    }

    /**
     * Based on the actual sampling. Not for MIS.
     */
    Rgb partial_light_flux(int node_index) {
        if (num_nodes <= 1) return DefaultRgb::black;
        if (node_index < 0 || node_index >= num_nodes) return DefaultRgb::black;

        return nodes[node_index].flux;
    }

    Rgb partial_eye_flux(int node_index) {
        if (num_nodes <= 1) return DefaultRgb::black;
        if (node_index < 0 || node_index >= num_nodes) return DefaultRgb::black;
        
        return nodes[node_index].importance;
    }

    Rgb radiance() {
        // we expect partial_light_flux(end) equal to partial_eye_flux(start) * Le. 
        // TODO: now it is not exactly equivalent because the cosine at eye is not considered properly (the last form factor when computed as light flux and eye flux)

        if (! is_complete) return DefaultRgb::black;

        if (is_unidirectional()) {
            if (nodes[0].is_on_light())
                return partial_light_flux(num_nodes - 1);
            else
                return partial_eye_flux(num_nodes - 1);
        } else {
            PathVertex &light_vertex = nodes[connection_index];
            PathVertex &eye_vertex = nodes[connection_index + 1];
            
            Rgb light_brdf = DefaultRgb::white;
            if (light_vertex.is_on_surface())
                light_brdf = light_vertex.material->get_bsdf()->eval(light_vertex.dg.uvn, light_vertex.wo, light_vertex.wi);
            Rgb eye_brdf = DefaultRgb::white;
            if (eye_vertex.is_on_surface())
                eye_brdf = eye_vertex.material->get_bsdf()->eval(eye_vertex.dg.uvn, eye_vertex.wo, eye_vertex.wi);
            Float G = FormFactor::form_factor_abs(light_vertex.dg.p, light_vertex.dg.n, eye_vertex.dg.p, eye_vertex.dg.n);
            
            Rgb radiance = 
                partial_light_flux(connection_index) *
                partial_eye_flux(connection_index + 1) * G * light_brdf * eye_brdf;

            return radiance;
        }
    }

    Rgb throughput() {

        if (num_nodes <= 1) return DefaultRgb::black;
        return nodes[num_nodes - 1].throughput;

    }

    /**
     * Probablity of how the path is generated 
     */
    Float pdf() {
        Float pdf_dA = 1.0f;
        for (int i = 0; i < num_nodes; ++i) {     // this already accounts for bidirectional path
            pdf_dA *= nodes[i].pdf_dA;
        }
        return pdf_dA;
    }

    static bool join(Scene *scene, 
                     MutablePath &light_subpath, MutablePath &eye_subpath, 
                     MutablePath &path) {

        Aggregate *agg = scene->get_aggregate();
        Float tmin = scene->get_tmin();
        Float max_tmax = scene->get_max_tmax();
        Float tick = scene->get_tick();

        HitRecord rec;
        PathVertex *light_vertex = &light_subpath.nodes[light_subpath.num_nodes - 1];
        PathVertex *eye_vertex = &eye_subpath.nodes[eye_subpath.num_nodes - 1];

        if (agg->hit(light_vertex->dg.p, eye_vertex->dg.p, tmin, tick)) {
            path.is_valid = false;
            return false;
        }

        path = light_subpath;

        int count = light_subpath.num_nodes;
        for (int i = eye_subpath.num_nodes - 1; i >= 0; --i) {
            path.nodes[count++] = eye_subpath.nodes[i];
        }

        path.connection_index = light_subpath.num_nodes - 1;

        light_vertex = &path.nodes[path.connection_index];
        eye_vertex = &path.nodes[path.connection_index + 1];

        light_vertex->wi = unit_vector(eye_vertex->dg.p - light_vertex->dg.p);
        eye_vertex->wi = -light_vertex->wi;

        light_vertex->pdf_by_next = pdf_area(*eye_vertex, *light_vertex);
        eye_vertex->pdf_by_next = pdf_area(*light_vertex, *eye_vertex);

        // pdf_by_prev is the actual pdf
        // pdf_by_next is the hypothesis pdf

        return true;
    }

};

typedef vector<MutablePath> MutablePaths;

/**
 * Iterate path starting from the first eye vertex
 */
class EyePathIterator {
public:
    EyePathIterator(const MutablePath &path) : path(path), reverse(false), index(-1) {
        if (path.num_nodes > 0 && path.nodes[0].is_on_light()) {
            reverse = true;
            index = path.num_nodes;
        }
    }

    bool has_next() const {
        if (! reverse) {
            if (index < path.num_nodes - 1) return true;
            else return false;
        }
        else {
            if (index > 0) return true;
            else return false;
        }
    }

    const PathVertex *current() {
        if (index < 0 || index >= path.num_nodes) return NULL;
        return &path.nodes[index];
    }

    void next() {
        if (path.num_nodes <= 0) return;

        if (! reverse) index++;
        else index--;
    }

private:
    const MutablePath &path;
    bool reverse;                   // when a path is started from light, reverse traversal is needed. 
    int index;
};

}    // end namespace 

#endif   