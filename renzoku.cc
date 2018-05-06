#include <GL/glew.h>
#include "common.h"

#include "scene.h"
#include "image.h"

#include "camera.h"
#include "shape.h"
#include "quad.h"
#include "triangle.h"
#include "sphere.h"
#include "surface.h"

#include "aggregate.h"

#include "lambertian.h"
#include "glass.h"
#include "mirror.h"
#include "phong.h"
#include "ward.h"
#include "thin_transparent.h"

#include "area_light.h"
#include "point_light.h"
#include "dir_light.h"

#include "sampler.h"
#include "stratified_sampler.h"

#include "first_bounce.h"
#include "monte_carlo_integrator.h"
#include "vpl.h"

#include "vpl_generator.h"

#include "obj_importer.h"
#include "json_importer.h"
#include "render_config_importer.h"

#include "mesh_view.h"
#include "image_block_view.h"
#include "viewer.h"
#include "glut_viewer.h"
//#include "qt_viewer.h"

#include "file.h"

using namespace Renzoku;
using Renzoku::Vec2;
using Renzoku::Vec3;

#include <iostream>
#include <vector>
#include <float.h>
using namespace std;

static void render_progressive(Scene *scene, int argc, char **argv) {    
    string name = scene->get_name();
    Size2 img_size = scene->get_image_size();
    int height = img_size.height;
    int width = img_size.width;
    
    Viewer *vr = new GlutViewer(height, width);
    
    /*
    QApplication *app = new QApplication(argc, argv);

    QGLFormat fmt;
    fmt.setDoubleBuffer(true);
    fmt.setRgba(true);
    fmt.setDepth(true);
    // my OpenGL code mixes between new OpenGL 3 and old OpenGL 2. Therefore, 
    // a compatibility profile of OpenGL 3.2 is needed. Sadly, from OS X 10.7, it seems
    // only OpenGL Core profile 3.2 is available. The legacy supported profile is 2.1.
    //fmt.setVersion(3, 2);    
    //fmt.setProfile( QGLFormat::CompatibilityProfile ); 
    //fmt.setProfile(QGLFormat::CoreProfile);
    fmt.setSampleBuffers(false);

    Viewer *vr = new QtViewer(app, fmt, height, width);
    ((QtViewer *)vr)->makeCurrent();      // make GL context before init GLEW
    */

    vr->set_scene(scene);     
    vr->init(argc, argv);
    
    
    /** 
     * GLEW needed for OpenGL shaders, FBO, etc.
     */    
    glewExperimental = GL_TRUE; // must be turned on for framebuffers EXT in core profile to work
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        string str((char *)glewGetErrorString(err));
        Log::critical(str, ExitCode::EXIT_CODE_LIBRARY_FAILED_TO_LOAD);
    }
    Log::info() << "GLEW: " << glewGetString(GLEW_VERSION) << endn;
               
    ParamDict *params = scene->get_param_dict();
    params->add("num_threads", 1);
    
    const int max_bounce = 2;       // set to 2 for Metropolis sampler
                                    
    // Guided path tracing with clustered VPL
    
    VirtualPointLight *integrator = new VirtualPointLight();
    scene->set_integrator(integrator);
    integrator->set_max_bounce(max_bounce);
    integrator->set_num_samples(10000);
    integrator->set_direct_lighting(false);
    integrator->set_radiance_estimate(VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE);
    integrator->set_start_particles(params->get_int("start_particles", 10000));
    integrator->set_num_passes(1);
    //integrator->set_radiance_estimate(VplRadianceEstimate::NONE);
    //integrator->set_radiance_estimate(VplRadianceEstimate::MRCS);
    //integrator->mrcs.set_max_clusters(300);       
    integrator->lightslice.set_num_clusters(params->get_int("clusters", 100));
    integrator->lightslice.set_max_clusters_per_slice(2 * params->get_int("clusters", 100));
    integrator->lightslice.set_max_slice_size(params->get_int("slice_size", 256));
    integrator->lightslice.set_density_radius(params->get_float("density_radius", 0.1));
    integrator->lightslice.set_num_neighbor_slices(32);
    integrator->set_virtual_indirect_samples(1);
    

	
    // Guided path tracing with Vorba's Gaussian fit
    /*
    VirtualPointLight *integrator = new VirtualPointLight();
    scene->set_integrator(integrator);
    integrator->set_max_bounce(max_bounce);
    integrator->set_num_samples(10000);
    integrator->set_direct_lighting(false);
    integrator->set_radiance_estimate(VplRadianceEstimate::VORBA_GUIDED_PATH_TRACE);
    integrator->set_virtual_indirect_samples(1);
    integrator->set_start_particles(params->get_int("start_particles", 10000));
    integrator->set_num_passes(1);
    integrator->set_save_every(params->get_int("save_every"));
    integrator->set_max_vpls(params->get_int("max_vpls"));
    */
	
    // Photon mapping
    /*
	VplPpm *integrator = new VplPpm();
    scene->set_integrator(integrator);
    integrator->set_max_bounce(max_bounce);
    integrator->set_num_samples(10000);
    integrator->set_direct_lighting(false);
	integrator->set_start_particles(100000);
	integrator->set_num_passes(100);
	integrator->set_save_every(params->get_int("save_every"));
	integrator->set_max_vpls(params->get_int("max_vpls"));

	integrator->use_progressive_photon(true);
	integrator->set_nearest_photons(64);
	integrator->set_alpha(0.67f);

    */

    /*
    NOTE: two sided surfaces is now misused to ensure smooth normals and wo is in the same direction
    for glass material to work. 
    This should be separated in fact. 
    */

    /*
    PathTracing *integrator = new PathTracing();
    //BidirPathTracing *integrator = new BidirPathTracing();
    //VertexConnectionMerging *integrator = new VertexConnectionMerging();    
    //BundlePathTracing *integrator = new BundlePathTracing();        
    //FirstBounce *integrator = new FirstBounce();
    //AmbientOcclusion *integrator = new AmbientOcclusion();
    //integrator->set_radius_scale(0.5f);
    //integrator->set_ray_samples(1);
    //LightTracing *integrator = new LightTracing();
    //PathPairing *integrator = new PathPairing();
    scene->set_integrator(integrator);    
    integrator->set_max_bounce(max_bounce);
    integrator->set_num_samples(10000);
    integrator->set_direct_lighting(false); 
    //integrator->set_share_light_path(true);
    //integrator->set_radius(16.0f); // small radius for speed
    scene->set_output_folder("data/path");
    */
    
    Log::info() << "Integrator: " << scene->get_integrator()->get_suffix() << endn;
    scene->print_statistics();
    
    scene->set_viewer(vr);
    if (scene->get_integrator()->is_viewer_outputing()) {
        Log::error() << "GPU view not supported." << endn;
        /*
        GpuView *gpu_view = new GpuView(scene, Size2(width, height));
        vr->add(gpu_view);
        */
    } else {
        ImageBlockView *image = new ImageBlockView(scene, Size2(width, height));
        scene->get_frame_buffer()->observable.attach_observer(image);
        vr->add(image);
    }
    MeshView *mesh = new MeshView(scene, height, width);    
    vr->add(mesh);

    BoundingBoxes boxes;
    scene->get_aggregate()->get_bounding_boxes(boxes);
    mesh->set_bounding_boxes(boxes);

    vr->show(0);
    vr->run();
}

void usage() {
    cout << "Usage:" << endl;
    cout << "renzoku <scene>" << endl;
}

int main(int argc, char **argv) {

    // Visual Studio: break if NaN occurs
    // this can cause performance penalty
    
    // turn on _EM_OVERFLOW might result in exceptions in half.h by OpenEXR 
    // because half is a 16-bit floating point representation which can only 
    // represents numbers between 6.1e-5 and 6.5e+4. See half.h.
    /*
    _clearfp();
    _controlfp(_controlfp(0, 0) & ~(_EM_INVALID | _EM_ZERODIVIDE),
               _MCW_EM);
    */

    // log all messages    
    Log::instance()->set_display_time(true);

	if (argc <= 1) {
		usage();
		return 0;
	}

	// make log file name
    std::time_t now = std::time(NULL);
    std::tm * ptm = std::localtime(&now);
    char buffer[32];        
    std::strftime(buffer, 32, "%m%d_%H%M%S", ptm); 
    ostringstream oss;
    string file_name, scene_path, scene_name;
    File::split_path(argv[1], file_name, scene_path, scene_name);
    oss << scene_name << "_" << buffer << ".txt";
    //Log::open(oss.str());

    if (argc == 2) {

        JsonImporter importer;
        importer.load(argv[1]);
        Scene *scene = importer.get_scene();
        //render_offline(scene);
        render_progressive(scene, argc, argv);    	
	    delete scene;

    } else if (argc == 3) {

        JsonImporter importer;
        importer.load(argv[1]);
        RenderConfigImporter integrator_importer;
        integrator_importer.load(argv[2]);

        Scene *scene = importer.get_scene();
        Integrator *integrator = integrator_importer.get_integrator();
        scene->set_integrator(integrator);

        render_progressive(scene, argc, argv);
        delete scene;

    } else {

        usage();

    }    
    return 0;
}
