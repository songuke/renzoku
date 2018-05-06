#ifndef _SCENE_H_
#define _SCENE_H_

#include "common.h"
#include "size.h"
#include "gl_vertex.h"
#include "pdf.h"
#include "param_dict.h"
#include "surface.h"
#include "aggregate.h"
#include "camera.h"

namespace Renzoku {

/**
 * This enumeration lists various ways to choose a light to sample 
 * when the scene has more than one light source. 
 * 
 * For diffuse scenes, power sampling can be a good strategy. 
 * 
 * In specular scenes, sometimes low-powered lights can cause specular highlights,
 * so they need to be sampled more. In such cases, uniform sampling is more preferred. 
 * For example, in Veach's MIS scene, uniform sampling is better.
 */
enum LightSelectionType {
    UNIFORM,
    POWER
};

struct SensorOptions {
    enum MultiSampleType {
        NONE,
        FSAA,
        MSAA
    };

    MultiSampleType multi_sample;
    int num_antialiased_samples;
};

struct ShadowMapOptions {
    Float polygon_offset_factor;
    Float polygon_offset_unit;
};

class Scene {
public:
    Scene(string name, 
          Aggregate *agg, Shapes *shapes, 
          Materials *materials, Textures *textures, 
          Surfaces *surfaces, 
          Lights *lights, Camera *camera, 
          Random *rd);

    ~Scene();
        
    inline Size2 get_image_size() const;

    inline BoundingBox get_bounding_box() const;
    inline BoundingSphere get_bounding_sphere() const;

    void set_light_selection_type(LightSelectionType type);

    /**
     * Return the index of the light that has a particular cdf.
     */
    int sample_light(Float r, Float &pdf);
    
    /**
     * Return the probability of sampling light k.
     */
    Float pdf_light(int k) const;
    Float pdf_light(Light *light) const;

    /**
     * Return a surface 
     */
    Surface *sample_surface(Float r, Float &pdf);
    
	/**
     * Return average albedo calculated by compute_average_albedo() function.
     * 
     * Can be used as a threshold to continue/stop in Russian roulette.
     */
	Float get_average_albedo() const;	

protected:
    /**
     * Generate a sampling distribution using power of lights in the scene. 
     *
     * The distribution can be used to select light with probability proportional to 
     * its emitted power. 
     */
    void generate_light_selection_pdf(LightSelectionType type);

	/**
     * Calculate the average value of all diffuse components of the materials in the scene.
     */
	void compute_average_albedo();
    
    /**
     * Generate a sampling distribution based on area to pick surfaces.
     */
    void generate_surface_selection_pdf();

public:
    inline Aggregate*  get_aggregate() const;
    inline Lights*     get_lights() const;
    inline Random*     get_random() const;
    inline Camera*     get_camera() const;
    inline Surfaces*   get_surfaces() const;
    inline Materials*  get_materials() const;
    inline Textures*   get_textures() const;

    inline void        set_tmin(Float tmin);
    inline Float       get_tmin() const;
    inline Float       get_max_tmax() const;
    inline Float       get_tick() const;
    
    inline string      get_name() const;    
    inline void        set_name(string s);
    string             get_auto_filename(const string &suffix = "", bool with_time_suffix = false) const;
    string             get_output_folder() const;
    void               set_output_folder(const string &path);

    void               set_camera(Camera *c);

    /**
     * Add a surface to a pre-loaded scene
     *
     * For adding parametric objects into scenes loading from OBJ files.
     */
    void        add_surface(Surface *surface);

    /**
     * Add a light to a pre-loaded scene.
     *
     * NOTE: area light source is a surface and must be added using add_surface() function
     * so that ray intersection test can be executed on area light sources. 
     */
    void        add_light(Light *light);

    /** 
     * Add an environment light.
     */
    void        add_env_light(ImageFloat *envmap);
        
    /**
     * Return all vertex coordinates of all surface triangles. 
     * Vertices of area light sources are excluded. 
     */
    void        get_triangle_coordinates(vector<Float> &coords);

    /**
     * Export all information for all vertices which can be used for OpenGL rendering. 
     */
    void        get_vertex_data(vector<GlVertex> &vertices, vector<int> &chunk_offsets);

    /**
     * Return the number of triangles, excluding the triangles from area light sources.
     */ 
    int         count_triangles();
    
    void        print_statistics();

    /**
     * Save geometry into a PLY file.
     */
    void        save_geometry(const string ply);
    void        save_geometry(const string ply, vector<Vec3> more_points);

    /** 
     * To control multisampling in GPU implementation.
     */
    inline SensorOptions     get_sensor_options() const;
    inline void              set_sensor_options(const SensorOptions &options);
    /**
     * Polygon offset and unit are scene dependent, hence make it a property of the scene.
     */
    inline ShadowMapOptions  get_shadow_map_options() const;
    inline void              set_shadow_map_options(const ShadowMapOptions &options);

    void                     set_integrator(Integrator *integrator);
    inline Integrator*       get_integrator() const;

    inline FrameBuffer*      get_frame_buffer() const;

    inline Stats*            get_stats_counter() const;

    inline void              use_facet_normals(bool facet);

           void              set_two_sided_surfaces(bool two);
    inline bool              get_two_sided_surfaces() const;

    /**
     * Only call this function if vertex normals are not available in the scene file.
     */
    void compute_vertex_normals();

    inline ThreadController* get_thread_controller() const;
    inline ParamDict*        get_param_dict() const;

    inline void         set_mesh_view(MeshView *mesh_view);
    inline MeshView*    get_mesh_view() const;

    inline void            set_image_view(ImageBlockView *image_view);
    inline ImageBlockView* get_image_view() const;

    inline void         set_viewer(Viewer *viewer);
    inline Viewer*      get_viewer() const;

    void start();

protected:
    Shapes*       shapes;
    Materials*    materials;
    Textures*     textures;
    Surfaces*     surfaces;                 // a general collection of surfaces, including area light sources. 
    Aggregate*    agg;
    Lights*       lights;                   // a general collection of all lights in the scene, including area lights, point lights, and environment lights.
    Camera*       camera;
    Random*       rd;
    Integrator*   integrator;
    FrameBuffer*  frame;                    // store the output of the integration

    // ray self-intersection constant
    Float tmin;
    Float max_tmax; 
    
    // time
    Float tick;

	Float         material_avg_albedo;
    
    // light selection
    LightSelectionType  light_selection_type;
    DiscretePdf         light_pdf;
    
    // surface selection
    DiscretePdf         surface_pdf;    

    string name;
    string output_folder;

    int num_triangles;    
    
    SensorOptions       sensor_options;
    ShadowMapOptions    shadow_map_options;

    Stats *stats_counter;
    
    bool facet_normals;
    bool two_sided;

    ThreadController *thread_controller;
    ParamDict *param_dict;
    
    Viewer *viewer;                         // for controlling and switching between views
    MeshView *mesh_view;                    // for visualizing debugging info such as bounding boxes
    ImageBlockView *image_view;             // for catching mouse/keyboard events in integrators
};

inline Random* Scene::get_random() const { 
    return rd; 
}

inline Aggregate* Scene::get_aggregate() const { 
    return agg; 
}

inline Lights* Scene::get_lights() const { 
    return lights;
}

inline Camera* Scene::get_camera() const { 
    return camera;
}

inline Surfaces *Scene::get_surfaces() const {
    return surfaces;
}

inline Materials *Scene::get_materials() const {
    return materials;
}

inline Textures *Scene::get_textures() const {
    return textures;
}

inline Float Scene::get_tmin() const { 
    return tmin;
}

inline void Scene::set_tmin(Float tmin) {
    this->tmin = tmin;
}

inline Float Scene::get_max_tmax() const { 
    return max_tmax; 
}

inline Float Scene::get_tick() const {
    return tick;
}

inline string Scene::get_name() const {
    return name;
}

inline void Scene::set_name(string name) {
    this->name = name;
}

inline string Scene::get_output_folder() const {
    return output_folder;
}

inline void Scene::set_output_folder(const string &path) {
    this->output_folder = path;
}

inline SensorOptions Scene::get_sensor_options() const {
    return sensor_options;
}

inline void Scene::set_sensor_options(const SensorOptions &options) {
    sensor_options = options;
}

inline Integrator* Scene::get_integrator() const {
    return integrator;
}

inline FrameBuffer* Scene::get_frame_buffer() const {
    return frame;
}

inline Stats* Scene::get_stats_counter() const {
    return stats_counter;
}

inline void Scene::use_facet_normals(bool facet) {
    facet_normals = facet;
}

inline ThreadController *Scene::get_thread_controller() const {
    return thread_controller;
}

inline ShadowMapOptions Scene::get_shadow_map_options() const {
    return shadow_map_options;
}

inline void Scene::set_shadow_map_options(const ShadowMapOptions &options) {
    shadow_map_options = options;
}

inline ParamDict* Scene::get_param_dict() const {
    return param_dict;
}

inline bool Scene::get_two_sided_surfaces() const {
    return two_sided;
}

inline void Scene::set_mesh_view(MeshView *mesh_view) {
    this->mesh_view = mesh_view;
}

inline MeshView *Scene::get_mesh_view() const {
    return mesh_view;
}

inline void Scene::set_image_view(ImageBlockView *image_view) {
    this->image_view = image_view;
}

inline ImageBlockView *Scene::get_image_view() const {
    return image_view;
}

inline void Scene::set_viewer(Viewer *viewer) {
    this->viewer = viewer;
}

inline Viewer *Scene::get_viewer() const {
    return viewer;
}

inline Size2 Scene::get_image_size() const {
    return camera->get_image_size();
}

/**
 * NOTE: currently each shape, e.g., triangle, stores its own bounding box and sphere. 
 * If memory is limited, we can change to compute bounding box/sphere on the fly.
 */
inline BoundingBox Scene::get_bounding_box() const {
    return agg->get_bounding_box();
}

inline BoundingSphere Scene::get_bounding_sphere() const {
    return agg->get_bounding_sphere();
}

};

#endif
