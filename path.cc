#include "path.h"

#include "random.h"
#include "material.h"

#include "light.h"
#include "area_light.h"
#include "point_light.h"
#include "spot_light.h"
#include "dir_light.h"
#include "concentric_mapping.h"

#include "shape.h"
#include "surface.h"
#include "aggregate.h"

#include "scene.h"
#include "shading_geometry.h"

namespace Renzoku {

PathNodePool *PathNodePool::_pool = NULL;

void SubpathGenerator::trace_eye_path(Scene *scene, const Ray &ray, int eye_verts, Path &path) {
    path.num_nodes = 0;
    if (eye_verts <= 0) {
        path.is_valid = false;
        return;
    }

    Rgb contrib = DefaultRgb::white;
    Float acc_pdf_dA = 1.0f;
    
    Ray r = ray;

    PathNode node(r.org(), r.dir(), DefaultRgb::white, 1.0f);
    path.set_last_node(node);

    Lights &lights = *scene->get_lights();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    
	// vertex 0 is the eye, so start with vertex 1    
    for (int vertex = 1; vertex < eye_verts; ++vertex) {        
        Float tmin     = scene->get_tmin();
        Float max_tmax = scene->get_max_tmax();
        Float tick     = scene->get_tick();
        
        HitRecord rec;
        if (! agg->hit(r, tmin, max_tmax, tick, rec)) break;
    
        LocalGeometry dg(r, rec);
        ShadingGeometry sg(r, rec);

        // second vertex considered sampled from eye
        acc_pdf_dA *= fabs(dot(dg.n, r.dir())) / (rec.t * rec.t);
        
        if (rec.light) {
            // traditional path tracing                      
            Rgb Le;
            rec.light->query(dg.p, dg.n, -r.dir(), Le);
            contrib *= Le;

		    // mark the last node of the path is a light
		    PathNode new_node(PathNode::EYE, dg.p, dg.n, rec.light, contrib, acc_pdf_dA);
		    path.set_last_node(new_node);    
            path.is_complete = true;
		    break;
        }
    
        Material *m = rec.material;                

        PathNode new_node(PathNode::EYE, dg, m, contrib, acc_pdf_dA);
	    path.set_last_node(new_node);
            
        Vec3 wi;
        Float pdf_wi;

        Vec3 wo = -r.dir();
        Rgb brdf = m->sample(rd, dg, wo, wi, pdf_wi);
        if (pdf_wi <= 0.0f) break;

        // TRICK: the local sampling already accounts for the angle far away from the out-going ray, 
        //        so the cosine here needs to account for the angle at the sampling origin.
        contrib     *= brdf * fabs(dot(dg.n, wi)) / pdf_wi;
        acc_pdf_dA  *= pdf_wi;

        r = Ray(dg.p, wi);
        node = new_node;
    }
}

void SubpathGenerator::trace_light_path(Scene *scene, int num_light_nodes, Path &path) {
    path.num_nodes = 0;
    if (num_light_nodes <= 0) {
        path.is_valid = false;
        return;
    }
    
    Random &rd = *scene->get_random();
    Vec2 sample; sample.random(rd);
    Float light_pdf;
    int k = scene->sample_light(sample.x(), light_pdf);
    if (k < 0) {
        path.is_valid = false;
        return;
    }

    Lights &lights = *scene->get_lights();
    Light *emitter = lights[k];
    
    Rgb contrib = DefaultRgb::white;
    Float acc_pdf_dA = light_pdf;

    if (emitter->get_light_type() == Light::AREA_LIGHT) {
        AreaLight *area_light = (AreaLight *)emitter;
        
        Vec3 q, qn, wo_q;
        Float pdf_q, pdf_q_wo;
        Rgb Le;                
        area_light->sample(rd, q, qn, wo_q, Le, pdf_q, pdf_q_wo);
    
        contrib     = DefaultRgb::white / pdf_q / light_pdf;
        acc_pdf_dA *= pdf_q;

        PathNode node(PathNode::LIGHT, q, qn, emitter, contrib, acc_pdf_dA);
        path.set_last_node(node);
        		
        Lights &lights = *scene->get_lights();
        Aggregate *agg = scene->get_aggregate();
        Random &rd     = *scene->get_random();

        Ray r(q, wo_q);
        contrib      = Le * fabs(dot(qn, wo_q)) / pdf_q / pdf_q_wo / light_pdf;        // FIXME: add light_pdf here to contrib
        acc_pdf_dA  *= pdf_q_wo;

        for (int vertex = 1; vertex < num_light_nodes; ++vertex) {
            Float tmin     = scene->get_tmin();
            Float max_tmax = scene->get_max_tmax();
            Float tick     = scene->get_tick();
        
            HitRecord rec;
            if (! agg->hit(r, tmin, max_tmax, tick, rec)) break;
               
            if (rec.light) break;               // consider light subpath ends at previous vertex            
            LocalGeometry dg(r, rec);

            Material *m = rec.material;
               	    
            acc_pdf_dA *= fabs(dot(dg.n, r.dir())) / (rec.t * rec.t);

            PathNode new_node(PathNode::LIGHT, dg, m, contrib, acc_pdf_dA);
	        path.set_last_node(new_node);
    
            Vec3 wi;
            Float pdf_wi;
            Rgb brdf = m->sample(rd, dg, -r.dir(), wi, pdf_wi);
            if (pdf_wi <= 0.0f) break;

            // TRICK: the local sampling already accounts for the angle far away from the out-going ray, 
            //        so the cosine here needs to account for the angle at the sampling origin.
            contrib     *= brdf * fabs(dot(dg.n, wi)) / pdf_wi;            
            acc_pdf_dA  *= pdf_wi;

            r = Ray(dg.p, wi);
            node = new_node;
        }
    }
}

bool SubpathGenerator::sample_direct_vertex(Scene *scene, PathNode &node, 
                                            DirectVertexSamplingRecord &sr) {
    node.acc_pdf_dA = 0.0f;
    sr.vertex1_shape = NULL;
    sr.vertex0_pdf_dA = 0.0f;
    sr.vertex1_pdf_dA = 0.0f;

    Random &rd = *scene->get_random();
    Vec2 sample; sample.random(rd);
    Float light_pdf;
    int k = scene->sample_light(sample.x(), light_pdf);
    if (k < 0 || k >= scene->get_lights()->size()) {        
        return false;
    }

    Lights &lights = *scene->get_lights();
    sr.emitter = lights[k];
    
    sr.vertex1_throughput = DefaultRgb::white;
    Rgb contrib = DefaultRgb::white;

    sr.vertex0_pdf_dA = light_pdf;
    Float acc_pdf_dA = light_pdf;

    Ray r;

    switch (sr.emitter->get_light_type()) {
        case Light::AREA_LIGHT: 
        {
            AreaLight *area_light = (AreaLight *)sr.emitter;
        
            Vec3 q, qn, wo_q;
            Float pdf_q, pdf_q_wo;
            Rgb Le;                
            area_light->sample(rd, q, qn, wo_q, Le, pdf_q, pdf_q_wo);
            if (pdf_q <= 0.0f || pdf_q_wo <= 0.0f) return false;

            acc_pdf_dA             *= pdf_q * pdf_q_wo;            
            contrib                 = Le * fabs(dot(qn, wo_q)) / acc_pdf_dA;
            sr.vertex1_throughput   = Le * fabs(dot(qn, wo_q));
            sr.vertex0_pdf_dA      *= pdf_q;
            sr.vertex1_pdf_dA       = pdf_q_wo;         // convert later

            Surface *surf = area_light->get_surface();
            sr.emitter_dg = LocalGeometry(q, qn, qn, surf->get_basis(q), 
                                          surf->get_shape()->tangent(), 
                                          surf->get_shape()->texture_uv(q)); 
            sr.emitter_dg_valid = true;
            sr.emitter_power = Le / sr.vertex0_pdf_dA;            

            r = Ray(q, wo_q);
            break;
        }

        case Light::POINT_LIGHT: 
        {
            PointLight *point_light = (PointLight *)(sr.emitter);
            Vec3 q;     
	        Vec3 wo;
	        Rgb Le;         
	        Float pdf_q, pdf_wo;
            point_light->sample(rd, q, wo, Le, pdf_q, pdf_wo);
            if (pdf_q <= 0.0f || pdf_wo <= 0.0f) return false;

            acc_pdf_dA  *= pdf_q * pdf_wo;
            contrib                 = Le / acc_pdf_dA;
            sr.vertex1_throughput   = Le;
            sr.vertex0_pdf_dA      *= pdf_q;
            sr.vertex1_pdf_dA       = pdf_wo;

            sr.emitter_power = Le / sr.vertex0_pdf_dA;
            sr.emitter_dg_valid = false;

            r = Ray(q, wo);
            break;
        }

        case Light::SPOT_LIGHT:
        {
            SpotLight *spot_light = (SpotLight *)(sr.emitter);
            Vec3 q;     
	        Vec3 wo;
	        Rgb Le;         
	        Float pdf_q, pdf_wo;
            spot_light->sample(rd, q, wo, Le, pdf_q, pdf_wo);
            if (pdf_q <= 0.0f || pdf_wo <= 0.0f) return false;

            acc_pdf_dA             *= pdf_q * pdf_wo;
            contrib                 = Le / acc_pdf_dA;
            sr.vertex1_throughput   = Le;
            sr.vertex0_pdf_dA      *= pdf_q;
            sr.vertex1_pdf_dA       = pdf_wo;

            sr.emitter_power = Le / sr.vertex0_pdf_dA;
            sr.emitter_dg_valid = false;

            r = Ray(q, wo);  
            break;
        }

        case Light::DIRECTIONAL_LIGHT: 
        {
            DirectionalLight *dir_light = (DirectionalLight *)(sr.emitter);
            Vec3 wi;
            Rgb Le;
            dir_light->query(wi, Le);

            // now generate VPLs to represent indirect illumination from directional light.            
            BoundingSphere bs = scene->get_aggregate()->get_bounding_sphere();
                       
            // uniform sampling on the disk using concentric mapping
            // then scale, rotate, translate.
            Vec3 proj_center = bs.center + 2.0f * bs.rad * wi;

            Vec2 d = square_to_disk(Vec2(rd(), rd())) * bs.rad;     
                    
            Vec3 n = wi;    // looking into negative z            
            Onb uvn;
            uvn.init_from_n(n);
            Vec3 p = uvn.local_to_world(Vec3(d.x(), d.y(), 0.0f)) + proj_center;
            
            acc_pdf_dA             *= 1.0f / (A_PI * bs.rad * bs.rad);
            contrib                 = Le / acc_pdf_dA;
            sr.vertex1_throughput   = Le;
            sr.vertex0_pdf_dA      *= 1.0f / (A_PI * bs.rad * bs.rad);
            sr.vertex1_pdf_dA       = 1.0f;

            sr.emitter_power = Le / sr.vertex0_pdf_dA;
            sr.emitter_dg_valid = false;

            r = Ray(p, -wi);
            break;
        }

        default: {
            return false;
        }
    }

    Aggregate *agg = scene->get_aggregate();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
        
    HitRecord rec;
    if (! agg->hit(r, tmin, max_tmax, tick, rec)) return false;
               
    if (rec.light) return false;               // consider light subpath ends at previous vertex            
    LocalGeometry dg(r, rec);

    Material *m = rec.material;
               	    
    // ensure T / acc_pdf_dA = contrib (the terms below actually cancels out)
    sr.vertex1_throughput *= fabs(dot(dg.n, r.dir())) / (rec.t * rec.t);
    sr.vertex1_pdf_dA *= fabs(dot(dg.n, r.dir())) / (rec.t * rec.t);

    acc_pdf_dA *= fabs(dot(dg.n, r.dir())) / (rec.t * rec.t);

    node = PathNode(PathNode::LIGHT, dg, m, contrib, acc_pdf_dA);
    node.actual_wi = -r.dir();

    sr.vertex1_shape = rec.shape;
    return true;
}

} // end namespace