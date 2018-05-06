#ifndef _PATH_H_
#define	_PATH_H_

#include "vec3.h"
#include "local_geometry.h"
#include "material.h"
#include "light.h"
#include "log.h"

namespace Renzoku {
    
struct Path;

struct PathNode {    
    LocalGeometry dg;

    Material *material;
    Light    *light;

    Rgb contrib[2];                 // equal to throughput / pdf that the path has so far                                    

    Float acc_pdf_dA;               // *accumulated* area measure
    
    int actual_prev;                // index parent node that generates this node    
    Vec3 actual_wi;    
    
    int path_index;                 // for retrieving the path that this node belongs to
    int index;

    enum Type {
        EYE,
        LIGHT
    } type;

    PathNode() : material(NULL), light(NULL), actual_prev(-1), type(EYE) {}

    /** 
     * Eye node
     */
    PathNode(const Vec3 &p, const Vec3 &wo_e, const Rgb &We, Float acc_pdf_dA) : 
        actual_prev(-1), 
        material(NULL), light(NULL),
        type(EYE), acc_pdf_dA(acc_pdf_dA) {
        
        contrib[EYE] = We;
        contrib[LIGHT] = DefaultRgb::black;
        dg.p = p;
        dg.n = wo_e;        // TODO: should this be the viewing vector instead? 
    }

    /**
     * General surface node.
     */
    PathNode(Type type, const LocalGeometry &dg, Material *m, const Rgb &T, Float acc_pdf_dA) :
        dg(dg), actual_prev(-1),
        material(m), light(NULL),
        type(type), acc_pdf_dA(acc_pdf_dA) {

        contrib[type] = T;
        contrib[1 - type] = DefaultRgb::black;
    }

    /**
     * Light node (starting a path or ending a path).
     */
    PathNode(Type type, const Vec3 &p, const Vec3 &n, Light *light, const Rgb &T, Float acc_pdf_dA) :
        actual_prev(-1),
        material(NULL), light(light), 
        type(type), acc_pdf_dA(acc_pdf_dA) {

        dg.p = p;
        dg.n = n;

        contrib[type] = T;
        contrib[1 - type] = DefaultRgb::black;
    }
    
    inline Rgb eye_contrib() const {
        return contrib[EYE];
    }

    inline Rgb light_contrib() const {
        return contrib[LIGHT];
    }

    inline bool is_on_light() const {
        return (light != NULL && material == NULL);
    }
    
    inline bool is_on_sensor() const {
        return (light == NULL && material == NULL && actual_prev == -1);
    }
};
typedef vector<PathNode> PathNodes;
typedef vector<PathNode *> PathNodePtrs;

const int MAX_PATH_NODE = 32;
const int DEFAULT_MAX_NODES = 262144 * MAX_PATH_NODE;
// TODO: work on path node pools for more coherency memory later. It is hard to allocate a big chunk of contiguous memory of this size.

class PathNodePool {
public:
    PathNodePool(int max_nodes = DEFAULT_MAX_NODES) : max_nodes(max_nodes) {
        nodes = new PathNode[max_nodes];

        for (int i = 0; i < DEFAULT_MAX_NODES; i += MAX_PATH_NODE) {
            free_list.push_back(nodes + i);
        }
    }

    ~PathNodePool() {
        if (nodes)
            delete[] nodes;
    }

    static PathNodePool *instance() {
        if (! _pool)
            _pool = new PathNodePool;
        return _pool;
    }
     
    PathNode *allocate(int size = MAX_PATH_NODE) {
        if (free_list.size() == 0) {
            // we should not relocate because several pointers are referring to the old buffer.
            Log::info() << "Buffer overflow." << endn;
            return NULL;
        }

        PathNode *nodes = free_list.back();
        free_list.pop_back();
        return nodes;
    }

    void clear(PathNode *nodes) {
        free_list.push_back(nodes);
    }

private:
    PathNode *nodes;    
    int max_nodes;

    vector<PathNode *> free_list;

    static PathNodePool *_pool;
};

struct Path {    
    PathNode nodes[MAX_PATH_NODE];
    int num_nodes;
    
    bool is_valid;                      /// when path is marked invalid, it is discarded.
    bool is_complete;                   /// indicates a full light path (with eye and light vertex found)

    /**
     * Create a path where its node buffer is provided from node pool.
     */
    //Path() : nodes(PathNodePool::instance()->allocate()), is_valid(true), is_complete(false), num_nodes(0) {
    Path() : is_valid(true), is_complete(false), num_nodes(0) {

    }

    ~Path() {
        //PathNodePool::instance()->clear(nodes);
    }

    Path(const Path &p) {
        memcpy(this->nodes, p.nodes, sizeof(PathNode) * MAX_PATH_NODE);
        this->num_nodes = p.num_nodes;
        this->is_valid = p.is_valid;
        this->is_complete = p.is_complete;
    }

    const Path *operator=(const Path &p) {
        memcpy(this->nodes, p.nodes, sizeof(PathNode) * MAX_PATH_NODE);
        this->num_nodes = p.num_nodes;
        this->is_valid = p.is_valid;
        this->is_complete = p.is_complete;
        return this;
    }
    
    PathNode &set_last_node(const PathNode &node) {
        this->nodes[num_nodes] = node;

        PathNode &last = this->nodes[num_nodes];
        last.actual_prev = num_nodes - 1;
        if (last.actual_prev >= 0)
            last.actual_wi = unit_vector(this->nodes[last.actual_prev].dg.p - last.dg.p);
        last.index = num_nodes;

        num_nodes++;
        return last;
    }

    const PathNode *get_last_node() const {
        if (num_nodes > 0)
            return &this->nodes[num_nodes - 1];
        return NULL;
    }
};
typedef vector<Path> Paths;

struct ExtendedPath {
    PathNode nodes[MAX_PATH_NODE];      /// last node of the path
    int num_nodes;
    bool is_valid;                      /// when path is marked invalid, it is discarded.
    bool is_complete;                   /// indicates a full light path (with eye and light vertex found)
    
    PathNode ancestors[2];              /// ancestor[0] is to generate the first node.
    
    ExtendedPath() : is_valid(true), is_complete(false), num_nodes(0) {

    }

    ExtendedPath(const ExtendedPath &p) {
        memcpy(this->nodes, p.nodes, sizeof(PathNode) * MAX_PATH_NODE);
        this->num_nodes = p.num_nodes;
        this->is_valid = p.is_valid;
        this->is_complete = p.is_complete;
    }

    const ExtendedPath *operator=(const ExtendedPath &p) {
        memcpy(this->nodes, p.nodes, sizeof(PathNode) * MAX_PATH_NODE);
        this->num_nodes = p.num_nodes;
        this->is_valid = p.is_valid;
        this->is_complete = p.is_complete;        
        return this;
    }
    
    void set_last_node(const PathNode &node) {
        this->nodes[num_nodes] = node;

        PathNode &last = this->nodes[num_nodes];
        last.actual_prev = num_nodes - 1;

        num_nodes++;
    }

    const PathNode *get_last_node() const {
        if (num_nodes > 0)
            return &this->nodes[num_nodes - 1];
        return NULL;
    }
};

struct DirectVertexSamplingRecord {
    // vertex on the emitter     
    Light *emitter;        
    Rgb emitter_power;
    LocalGeometry emitter_dg;
    bool emitter_dg_valid;

    // vertex on the surface illuminated by the emitter
    Shape *vertex1_shape;
    Rgb vertex1_throughput;

    Float vertex0_pdf_dA;               // area measure of vertex 0 and 1
    Float vertex1_pdf_dA;

    DirectVertexSamplingRecord() : 
        emitter_dg_valid(false), emitter(NULL), 
        vertex1_shape(NULL), vertex1_pdf_dA(0.0f), vertex0_pdf_dA(0.0f) {
    }
};

class SubpathGenerator {
public:
    static void trace_eye_path(Scene *scene, const Ray &ray, int eye_verts, Path &path);
    static void trace_light_path(Scene *scene, int num_light_nodes, Path &path);

    /**
     * Generate the path vertex that is directly visible from a light.
     *
     * This path vertex can be used as a VPL to render indirect illumination.
     */
    static bool sample_direct_vertex(Scene *scene, PathNode &node,
                                     DirectVertexSamplingRecord &sr);
};

} // end namespace
#endif

