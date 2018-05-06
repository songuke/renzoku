#include "aggregate.h"
#include "lambertian.h"

#include <stdexcept>

namespace Renzoku {
Aggregate::Aggregate(const Surfaces &_surfaces) : surfaces(_surfaces) {
    compute_bounding_box();
    compute_bounding_sphere();

    compute_density();    
    set_surface_id();
    //compute_average_albedo();
}

Aggregate::~Aggregate() {
    
}

Float Aggregate::area() const {
	Float sum = 0.;
	for(int k=0; k< surfaces.size(); k++) {
		sum += surfaces[k].area();
	}
	return sum;
}

BoundingBox Aggregate::get_bounding_box() const {
    return box;
}

BoundingSphere Aggregate::get_bounding_sphere() const {
    return sphere;
}

BoundingBox Aggregate::get_bounding_box(Surfaces &surfaces) {
    BoundingBox box;
    for (int i = 0; i < surfaces.size(); ++i) 
        box.merge(surfaces[i].get_bounding_box());
    return box;
}

void Aggregate::compute_bounding_box() {
    box.reset();
    for (int i = 0; i < surfaces.size(); ++i) 
        box.merge(surfaces[i].get_bounding_box());    
}

void Aggregate::compute_bounding_sphere() {
    sphere.reset();
    for (int i = 0; i < surfaces.size(); ++i) {
        sphere.merge(surfaces[i].get_bounding_sphere());
    }
}

bool Aggregate::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {    
    GeometryHit best_gh;
    for (int k = 0; k < surfaces.size(); ++k) {        
        GeometryHit gh;
        if (surfaces[k].get_shape()->hit(r, tmin, tmax, time, gh)) {
            tmax = gh.t;
            best_gh = gh;
        }
    }
    
    record.copy_geometry_hit(best_gh);
    if (best_gh.hit) {
        record.shape->fill_hit_record(r, best_gh, record);
        const Surface &s = surfaces[record.shape->get_surface_index()];
        s.fill_hit_record(r, record);
    }
    return best_gh.hit;
}

bool Aggregate::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    for (int k = 0; k < surfaces.size(); ++k) {
        if (surfaces[k].hit(r, tmin, tmax, time)) {
            return true;
        }
    }
    return false;
}

bool Aggregate::hit(const Vec3 &p, const Vec3 &q, Float tmin, Float time) const {
    Ray shadow(p, q - p);
    Float shadow_tmax = (q - p).length() - tmin;  // deduct an epsilon amount to avoid
	                                              // the case where q is 'behind' the surface of 
	                                              // itself.
    return this->hit(shadow, tmin, shadow_tmax, time);
}

void Aggregate::compute_density() {
    shape_pdf = new Float[surfaces.size()];
    shape_cdf = new Float[surfaces.size() + 1];

    shape_sum_area = 0.;
    for (int i = 0; i < surfaces.size(); ++i) {
        shape_sum_area += surfaces[i].area();
    }

    for (int i = 0; i < surfaces.size(); ++i) {
        if (surfaces[i].get_env_light() == NULL) { // don't sample an env light
            shape_pdf[i] = surfaces[i].area() / shape_sum_area;
        }
    }
    
    if (surfaces.size() > 0) {
        shape_cdf[0] = 0;
        for (int i = 1; i <= surfaces.size(); ++i) {
            shape_cdf[i] = shape_cdf[i - 1] + shape_pdf[i - 1];
        }
    }
}

void Aggregate::set_surface_id() {
    for (int i = 0; i < surfaces.size(); ++i) {
        surfaces[i].get_shape()->set_surface_index(i);
    }
}

int Aggregate::sample_surface(Random &rd, Float &pdf) const {
    Float u = rd();
    for (int i = 0; i < surfaces.size(); ++i) {
        if (u >= shape_cdf[i] && u < shape_cdf[i + 1]) {
            pdf = shape_pdf[i];
            return i;
        }
    }
    // cout << "Found 1." << endl;
    // u = 1 then return the final
    if (surfaces.size() > 0) {
        pdf = shape_pdf[(int)surfaces.size() - 1];
    }
    return surfaces.size() - 1;

}

void Aggregate::sample(Random &rd, Vec3 &p, Vec3 &n, Float &pdf) const {
    Float pdf_k;    
    int k = sample_surface(rd, pdf_k);
    if (k < 0) {
        // no surface exists
        throw std::runtime_error("No surfaces exist. Sampling failed.");

        // TODO: perhaps a better design is to return false. 
        return;
    }
    Shape *shape = surfaces[k].get_shape();
    Float pdf_shape;    
    shape->sample(rd, p, n, pdf_shape); 
    pdf = pdf_k * pdf_shape;
}

} // end namespace Renzoku

