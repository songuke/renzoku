#include "scene.h"

#include "material.h"
#include "textured_material.h"
#include "shape.h"
#include "surface.h"

#include "light.h"
#include "area_light.h"
#include "env_light.h"

#include "log.h"
#include "stats.h"

#include "sampler.h"
#include "integrator.h"
#include "integrator_scheduler.h"
#include "frame.h"
#include "thread_controller.h"

#include <sstream>
using namespace std;

namespace Renzoku {

Scene::Scene(string name,
             Aggregate *agg, 
             Shapes *shapes, Materials *materials, Textures *textures, Surfaces *surfaces, 
             Lights *lights, Camera *camera, 
             Random *rd) 
    : name(name), agg(agg), 
      shapes(shapes), materials(materials), textures(textures), 
      surfaces(surfaces), lights(lights), camera(camera),
      rd(rd)
{

    tmin        = RAY_TMIN;
    max_tmax    = RAY_TMAX;
    tick        = .0f;
    
    agg->set_scene(this);
    light_selection_type = POWER;
    this->generate_light_selection_pdf(light_selection_type);
        
	compute_average_albedo();
            
    for (int i = 0; i < materials->size(); ++i) {
        (*materials)[i]->set_material_index(i);
    }

    this->generate_surface_selection_pdf();

    this->sensor_options.multi_sample = SensorOptions::NONE;
    this->sensor_options.num_antialiased_samples = 4;  

    this->shadow_map_options.polygon_offset_factor  = 1.0f;
    this->shadow_map_options.polygon_offset_unit    = 4.0f;
    this->param_dict = new ParamDict();

    this->frame = new FrameBuffer(camera->get_image_size());
    this->stats_counter = new Stats;

    this->output_folder = ".";

    facet_normals = false;
    two_sided = false;
    thread_controller = new ThreadController();
    mesh_view = NULL;
}

Scene::~Scene() {
    for (Shapes::iterator i = shapes->begin(); i != shapes->end(); ++i) {
        Shape *s = *i;
        delete s;
    }
    for (Lights::iterator i = lights->begin(); i != lights->end(); ++i) {
        Light *l = *i;
        delete l;
    } 
    for (Materials::iterator i = materials->begin(); i != materials->end(); ++i) {
        Material *m = *i;
        delete m;
    }
    for (Textures::iterator i = textures->begin(); i != textures->end(); ++i) {
        delete (*i);
    }
    delete agg;
    delete camera;
    delete lights;
    delete materials;
    delete textures;
    delete shapes;
    delete surfaces;
    delete rd;
    delete frame;
    delete stats_counter;
}

void Scene::print_statistics() {    
    if (lights->size() <= 0) 
        Log::warn() << "No light source exists in the scene." << endn;
    else
        Log::info() << "Number of lights: " << lights->size() << endn;
}

void Scene::set_camera(Camera *c) {
    if (camera)
        delete camera;

    this->camera = c;
}

void Scene::set_integrator(Integrator *integrator) {
    this->integrator = integrator;
    integrator->observable.attach_observer(frame);
}

void Scene::add_surface(Surface *surface) {    
    surfaces->append(*surface);
    shapes->push_back(surface->get_shape());
    
    Material *m = surface->get_material();
    if (m) {
        materials->push_back(m);
    } else {
        AreaLight *area_light = surface->get_area_light();
        if (area_light) {
            this->add_light(area_light);
        } else {
            EnvLight *env_light = surface->get_env_light();
            if (env_light) {
                this->add_light(env_light);
            }
        }
    }
        
    // the list of surfaces should be maintained by scene. 
    // Aggregate should not modify this list.
    agg->add_surface(surface);
}

void Scene::add_light(Light *light) {
    lights->push_back(light);
    
    // refresh sampling distribution
    this->generate_light_selection_pdf(light_selection_type);
}

void Scene::add_env_light(ImageFloat *envmap) {
    BoundingBox bb = agg->get_bounding_box();

    if (camera == NULL) {
        Log::info() << "Environment light is modelled as a sphere and needs to bound the camera. Please add camera first." << endn;
        return;
    }

    // ensure camera is inside the bounding sphere
    BoundingBox camera_bb(camera->get_eye());
    bb.merge(camera_bb);

    Sphere *sphere = new Sphere(bb.centroid(), bb.size().max_component());

    Surface *env = new Surface(sphere, envmap);

    agg->add_surface(env);

    lights->push_back(env->get_env_light());
}

void Scene::set_light_selection_type(LightSelectionType type) {
    this->light_selection_type = type;

    // refresh sampling distribution
    this->generate_light_selection_pdf(type);
}

void Scene::generate_light_selection_pdf(LightSelectionType type) {
    vector<Float> dist;    
    dist.resize(lights->size());
    for (int i = 0; i < lights->size(); ++i) {
        switch (type) { 
        case POWER:
            dist[i] = (*lights)[i]->power().luminance();        
            break;

        case UNIFORM:            
            dist[i] = 1.0f;            
            break;
        }
    }
    light_pdf.set_distribution(dist);

    // store the sampling probability to each light for easy query in pdf_light(Light *light) function.
    for (int i = 0; i < lights->size(); ++i) {
        (*lights)[i]->set_selection_pdf(light_pdf.probability(i));
    }
}
    
int Scene::sample_light(Float rand, Float &pdf) {
    int k;
    light_pdf.sample(rand, k, pdf);
    return k;
}

Float Scene::pdf_light(int k) const {
    return light_pdf.probability(k);
}

Float Scene::pdf_light(Light *light) const {
    return light->selection_pdf();
}

void Scene::generate_surface_selection_pdf() {
    vector<Float> dist;    
    dist.reserve(surfaces->size());
    for (int i = 0; i < surfaces->size(); ++i) {
        if ((*surfaces)[i].get_material() != NULL) {
            dist.push_back((*surfaces)[i].area());
        }
    }
    surface_pdf.set_distribution(dist);
}

Surface *Scene::sample_surface(Float r, Float &pdf) {
    int k;
    surface_pdf.sample(r, k, pdf);
    if (k >= 0 && k < surfaces->size()) return &(*surfaces)[k];
    return NULL;
}

void Scene::compute_average_albedo() {
    if (surfaces->size() <= 0) return;

    material_avg_albedo = 0.;
	Float total_area = 0;
    for (int i = 0; i < surfaces->size(); ++i) {
        Surface &s = (*surfaces)[i];

        // use the diffuse component only, weight by area
        Material *m = s.get_material();
        if (m) {
            //Rgb albedo = m->get_diffuse_component();        
            Rgb albedo = m->get_representative_color();
            material_avg_albedo += albedo.avg() * s.area();
		    total_area += s.area();
        }
    }
    material_avg_albedo /= total_area;
}

Float Scene::get_average_albedo() const {
	return material_avg_albedo;
}

int Scene::count_triangles() {
    int num_triangles = 0;
    for (int i = 0; i < shapes->size(); ++i)
        num_triangles += (*shapes)[i]->get_triangle_count();
    return num_triangles;
}

void Scene::get_triangle_coordinates(vector<Float> &coords) {
    int num_triangles = count_triangles();

    coords.reserve(num_triangles * 3); // avoid reallocation
    for (int i = 0; i < shapes->size(); ++i) {
        (*shapes)[i]->get_triangle_coordinates(coords);
    }
}
   
void Scene::get_vertex_data(vector<GlVertex> &vertices, vector<int> &chunk_offsets) {
    int num_triangles = count_triangles();

    vertices.clear();
    chunk_offsets.clear();

    chunk_offsets.reserve(textures->size() + 1);
    vertices.reserve(num_triangles * 3);
        
    // NOTE: GLSL 4 still does not support dynamic index into an array of texture samplers,
    // so the safe approach is to render object one by one, each has its own set of textures.
    for (int j = 0; j < surfaces->size(); ++j) {
        Material *m = (*surfaces)[j].get_material();
        if (m == NULL || 
            m->get_material_type() != Material::TEXTURED) {

            (*surfaces)[j].get_vertex_data(vertices);
        }
    }
    chunk_offsets.push_back(vertices.size());
    for (int i = 0; i < textures->size(); ++i) {
        for (int j = 0; j < surfaces->size(); ++j) {
            Material *m = (*surfaces)[j].get_material();
            if (m != NULL &&
                m->get_material_type() == Material::TEXTURED &&
                ((TexturedMaterial *)m)->get_texture() == (*textures)[i]) {

                (*surfaces)[j].get_vertex_data(vertices);

            }
        }
        chunk_offsets.push_back(vertices.size());
    }
}

void Scene::save_geometry(const string ply) {
    vector<Vec3> v;
    save_geometry(ply, v);
}
    
void Scene::save_geometry(const string ply, vector<Vec3> more_points) {    
    // write surfaces: only triangle list is supported.    
    vector<GlVertex> vertices;
    vector<int> offsets;
    this->get_vertex_data(vertices, offsets);

    int num_triangles = count_triangles();

    ostringstream oss;
    oss << ply << ".ply";

    FILE *f = fopen(oss.str().c_str(), "w");
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "element vertex %d\n", vertices.size() + more_points.size());
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "element face %d\n", num_triangles);
    fprintf(f, "property list uchar int vertex_index\n");
    fprintf(f, "end_header\n");
    for (int i = 0; i < vertices.size(); ++i) {
        GlVertex &v = vertices[i];
        fprintf(f, "%f %f %f\n", v.pos[0], v.pos[1], v.pos[2]);
    }

    // more points
    for (int i = 0; i < more_points.size(); ++i) {
        Vec3 p = more_points[i];
        fprintf(f, "%f %f %f\n", p.x(), p.y(), p.z());
    }

    for (int i = 0; i < num_triangles; ++i) {
        fprintf(f, "%d %d %d %d\n", 3, 3*i, 3*i+1, 3*i+2);
    }
    
    fclose(f); 
}

struct Vertex {
    Vec3 p;
    int shape_idx;
    int vertex_idx;
};

static bool vertex_lessthan(const Vertex &v1, const Vertex &v2) {
    Float x1 = v1.p.x();
    Float y1 = v1.p.y();
    Float z1 = v1.p.z();
    Float x2 = v2.p.x();
    Float y2 = v2.p.y();
    Float z2 = v2.p.z();

	if (x1 < x2) 
        return true;
    else 
        if (x1 == x2) {
            if (y1 < y2)
                return true;
            else 
                if (y1 == y2) {
                    if (z1 < z2)
                        return true;     
                    else 
                        return false;
                } else return false;
        } else return false;
}

static bool vertex_equal(const Vertex &v1, const Vertex &v2) {
    return (v1.p.x() == v2.p.x() && v1.p.y() == v2.p.y() && v1.p.z() == v2.p.z());
}

void Scene::compute_vertex_normals() {
    vector<Vec3> positions;    
    positions.reserve(4);
    
    vector<Vertex> vertices;
    vertices.reserve(shapes->size() * 4);

    int k;
    Shapes::iterator i;
    for (k = 0, i = shapes->begin(); i != shapes->end(); ++i, ++k) {
        Shape *s = *i;
        positions.clear();        
        s->get_vertex_positions(positions);
        
        for (int j = 0; j < positions.size(); ++j) {
            Vertex v;
            v.p = positions[j]; 
            v.shape_idx = k;
            v.vertex_idx = j;
            vertices.push_back(v);
        }        
    }
    Log::info() << "Total vertices: " << vertices.size() << endn;
    std::sort(vertices.begin(), vertices.end(), vertex_lessthan);

    if (vertices.size() == 1) return;
    
    Vec3 normal = (*shapes)[vertices[0].shape_idx]->normal(vertices[0].p);
    int start = 0;
    int end = 1;
    while (end < vertices.size()) {
        while (vertex_equal(vertices[end], vertices[start])) {
            normal += (*shapes)[vertices[end].shape_idx]->normal(vertices[end].p);
            ++end;
        }
        normal = unit_vector(normal);

        for (int k = start; k < end; ++k) {
            (*shapes)[vertices[k].shape_idx]->set_vertex_normal(vertices[k].vertex_idx, normal);
        }

        if (end == vertices.size()) break;

        normal = (*shapes)[vertices[end].shape_idx]->normal(vertices[end].p);
        start = end;
        ++end;
    }
}

/**
 * Call by a view's on show.
 */
void Scene::start() {
    if (! integrator->is_viewer_outputing()) { // not a GPU integrator        
        thread_controller->kill_all();
    }

    integrator->initialize(this);

    if (! integrator->is_viewer_outputing()) {

        if (integrator->get_suffix() == "vmlt") {

            boost::thread *t = new boost::thread(boost::bind(&Integrator::integrate, integrator));
            thread_controller->add_thread(t);

        } else {

            IntegratorScheduler scheduler;
            scheduler.initialize(this);
        
            scheduler.set_num_threads(param_dict->get_int("num_threads", 1));

            boost::thread *t = new boost::thread(boost::bind(&IntegratorScheduler::start, scheduler));
            thread_controller->add_thread(t);

        }
        
    } else {
        // integrate_partial() is called by ImageView::display()
    }
}

string Scene::get_auto_filename(const string &suffix, bool with_time_suffix) const {
    
    ostringstream oss;
    oss << this->get_name() << "_" << integrator->get_suffix();
    
    if (suffix != "") 
        oss << "_" << suffix;

    if (with_time_suffix) {
        std::time_t now = std::time(NULL);
        std::tm * ptm = std::localtime(&now);
        char buffer[32];        
        std::strftime(buffer, 32, "%m%d_%H%M%S", ptm); 

        oss << "_" << buffer;
    }

    oss << ".exr";
    return oss.str();
}

void Scene::set_two_sided_surfaces(bool two) {
    two_sided = two;

    for (int i = 0; i < surfaces->size(); ++i) {
        (*surfaces)[i].set_two_sided(two);
    }
}

} // end namespace
