#include "render_config_importer.h"

#include "monte_carlo_integrator.h";
#include "path_tracing.h"
#include "vpl.h"
//#include "vpl_gpu.h"
//#include "vsl_gpu.h"
#include "lightslice.h"

namespace Renzoku {

string RenderConfigImporter::keywords[] = {
    "integrator"
};

RenderConfigImporter::RenderConfigImporter() {
}

bool RenderConfigImporter::process_object(int stack_level, string last_key, JsonObject *obj) {
    if (stack_level == 1) {
        if (last_key == keywords[KW_INTEGRATOR]) {
            Log::info() << "Load integrator: " << endn;
            load_integrator(obj);
            return true;
        }
        else
            return false;
    }
    return false;
}

void RenderConfigImporter::process_root(JsonObject *root) {
    if (root->is_empty()) {
        Log::info() << "Nothing left for further processing." << endn;
    }
}

bool RenderConfigImporter::process_array(string last_key, JsonArray *array) {
    return false;
}

Integrator *RenderConfigImporter::get_integrator() {
    return integrator;
}

void RenderConfigImporter::load(const string _file) {
    string config_name;
    string config_path;

    string file = _file;
    unsigned int found_dot = file.rfind(".json");
    if (found_dot == string::npos) {

        if (file[file.length() - 1] == '/') {
            file = file.substr(0, file.length() - 1);
        }
        unsigned int found_slash = file.find_last_of('/');
        config_name = file.substr(found_slash + 1);
        config_path = file + "/";
        file = config_path + config_name + ".json";

    }
    else {
        unsigned int found_slash = file.find_last_of('/');
        config_name = file.substr(found_slash + 1, found_dot - found_slash - 1);
        config_path = file.substr(0, found_slash + 1);
    }
    Log::info() << "Render config: " << file << endn;

    JsonParser parser(this);
    parser.load(file);
}

void RenderConfigImporter::load_integrator(JsonObject *spec) {
    int max_bounce = spec->get_value("max_bounce")->to_int();
    if (max_bounce == 0) {
        Log::warn() << "Integrator max bounce is zero." << endn;
    }

    string type = spec->get_value("type")->to_string();

    if (type == "path") {
        int samples = spec->get_value("samples")->to_int();
        if (samples == 0) {
            Log::warn() << "Monte Carlo integrator total sample is zero." << endn;
        }
        MonteCarloIntegrator *integrator = NULL;

        if (type == "path")
            integrator = new PathTracing();

        integrator->set_max_bounce(max_bounce);
        integrator->set_num_samples(samples);
        return;

    }

    else if (type == "vpl") {

        // move all into constructor VirtualPointLight::load_parameters(JsonObject *) after detecting the type of the integrator

        VirtualPointLight *integrator = new VirtualPointLight();
        integrator->set_max_bounce(max_bounce);
        integrator->set_direct_lighting(spec->get_bool("direct_lighting"));
        //integrator->set_num_samples(10000);

        VplRadianceEstimate::Type radiance_estimate_type;
        string radiance_estimate = spec->get_string("radiance_estimate");
        if (radiance_estimate == "lightslice_guided_path_trace") {
            radiance_estimate_type = VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE;
        } // more here
        integrator->set_radiance_estimate(radiance_estimate_type);

        switch (radiance_estimate_type) {
        case VplRadianceEstimate::LIGHTSLICE_GUIDED_PATH_TRACE:
            JsonObject *lightslice = spec->get_value_strict("lightslice")->to_object();
            integrator->lightslice.set_num_clusters(lightslice->get_int("clusters"));
            integrator->lightslice.set_max_clusters_per_slice(lightslice->get_int("max_clusters"));
            integrator->lightslice.set_max_slice_size(lightslice->get_int("slice_size"));
            integrator->lightslice.set_density_radius(lightslice->get_numeric("density_radius"));
            integrator->lightslice.set_num_neighbor_slices(lightslice->get_int("neighbor_slices"));
            break;
        }
            
        integrator->set_start_particles(spec->get_int("start_particles"));

        string clamp = spec->get_string("clamp_type");
        VplRadianceClamp::Type clamp_type;
        if (clamp == "radiance") {
            clamp_type = VplRadianceClamp::RADIANCE;
        } // more here
        integrator->set_clamping(clamp_type);

        JsonObject *clamp_threshold = spec->get_value_strict("clamp_threshold")->to_object();
        switch (clamp_type) {
        case VplRadianceClamp::RADIANCE:
            integrator->set_clamping_threshold(clamp_threshold->get_numeric("radiance"));
            integrator->set_incident_radiance_clamping_threshold(clamp_threshold->get_numeric("incident_radiance"));
            break;
        }
        integrator->set_clamping_relaxation(spec->get_numeric("clamp_relaxation"));

        this->integrator = integrator;
        return;
    }
    /*
    else if (type == "vpl_gpu") {
        int particles = spec->get_value("particles")->to_int();
        if (particles == 0) {
            Log::warn() << "Virtual point light integrator total starting particle is zero." << endn;
        }

        bool direct = spec->get_value("direct")->to_bool();

        VirtualSphericalLightGpu *integrator = new VirtualSphericalLightGpu();
        integrator->set_max_bounce(max_bounce);
        integrator->set_start_particles(particles);
        integrator->set_direct_lighting(direct);
        return;

    }
    else if (type == "vsl_gpu") {
        int particles = spec->get_value("particles")->to_int();
        if (particles == 0) {
            Log::warn() << "Virtual spherical light integrator total starting particle is zero." << endn;
        }

        bool direct = spec->get_value("direct")->to_bool();

        VirtualSphericalLightGpu *integrator = new VirtualSphericalLightGpu();
        integrator->set_max_bounce(max_bounce);
        integrator->set_start_particles(particles);
        integrator->set_direct_lighting(direct);
    }*/
    else {
        Log::error() << "Invalid integrator type: " << type << "." << endn;;
    }
}

}  // end namespace