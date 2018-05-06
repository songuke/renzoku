#include "light_transport_matrix.h"

#include "scene.h"
#include "shape.h"
#include "aggregate.h"
#include "camera.h"
#include "brdf_point_light.h"
#include "stats.h"

namespace Renzoku {

LightTransportMatrix::LightTransportMatrix(int rows, int cols) {
    m.resize(rows, cols);
}

void LightTransportMatrix::fill(Scene *scene, 
                                IVirtualPointLightEvaluator *evaluator, 
                                BrdfPointLights &vpls, 
                                Pixels &pixels) {

    Camera *camera = scene->get_camera();
    Aggregate *agg = scene->get_aggregate();
    Random &rd     = *scene->get_random();
    Float tmin     = scene->get_tmin();
    Float max_tmax = scene->get_max_tmax();
    Float tick     = scene->get_tick();
    Lights &lights = *scene->get_lights();
    
    m.resize(pixels.size(), vpls.size());

    for (int ridx = 0; ridx < pixels.size(); ++ridx) {
        Vec2 pixel = pixels[ridx];
        
        Ray r = camera->shoot_ray(pixel.y(), pixel.x());
        HitRecord rec;
        if (! agg->hit(r, tmin, max_tmax, tick, rec)) continue;

        // skip pixels that hit the light source directly since light source does not reflect light so 
        // VPLs won't affect the light source.
        if (rec.light) continue;

        Receiver recv(r, rec);
        for (int k = 0; k < vpls.size(); ++k) {                        
            Rgb radiance = evaluator->radiance(scene, recv, vpls[k]);
            (this->m)(ridx, k) = radiance;                 
        }
    }    
}

void LightTransportMatrix::get_column(int i, vector<Float> &c) const {    
    int rows = m.get_rows();
    if (c.size() != rows * 3)
        c.resize(rows * 3);

    for (int r = 0; r < rows; ++r) {
        Rgb v = m(r, i);
        c[3 * r    ] = v.red();
        c[3 * r + 1] = v.green();
        c[3 * r + 2] = v.blue();
    }
}

void LightTransportMatrix::set_column(int i, const vector<Float> &c) {
    int rows = m.get_rows();
    
    for (int r = 0; r < rows; ++r) {
        m(r, i) = Rgb(c[3 * r], c[3 * r + 1], c[3 * r + 2]);
    }
}


void LightTransportMatrix::generate_alphas_for_clustering(vector<Float> &alphas) {
    int rows = m.get_rows();
    int cols = m.get_cols();
    
    MatrixC<Float> R(rows, cols);
    m.flatten_magnitude(R);

    // alpha = (n_R^t * n_R - R^t * R) 1
    //       = (n_R^t * n_R * 1 - R^t * R * 1)

    MatrixC<Float> column_norms(1, cols);
    for (int j = 0; j < cols; ++j) {
        Float sum = 0.0f;
        for (int i = 0; i < rows; ++i) {
            sum += R(i, j) * R(i, j);            
        }
        column_norms(0, j) = sqrt(sum);        
    }
    
    MatrixC<Float> Rt(cols, rows);
    R.transpose(Rt);
        
    vector<Float> u(rows);
    R.row_sum(u);                       // u = R * 1
    vector<Float> v(cols);    
    Rt.mul(u, v);                       // v = Rt * u
    
    vector<Float> c;                    // n_R * 1
    column_norms.row_sum(c); 

    alphas.resize(cols);
    for (int i = 0; i < cols; ++i) {
        alphas[i] = column_norms(0, i) * c[0] - v[i];
    }
}

Float LightTransportMatrix::get_column_norm(int i) const {
    Float magnitude = 0.0f;
    for (int r = 0; r < m.get_rows(); ++r) {
        Rgb v = m(r, i);
        magnitude += v.red() * v.red() + v.green() * v.green() + v.blue() * v.blue();                        
    }
    return sqrt(magnitude);
}

void LightTransportMatrix::get_column_norms(vector<Float> &norms) const {
    norms.resize(m.get_cols());
        
    for (int c = 0; c < m.get_cols(); ++c) {
        Float magnitude = 0.f;

        for (int r = 0; r < m.get_rows(); ++r) {
            Rgb v = m(r, c);
            magnitude += v.red() * v.red() + v.green() * v.green() + v.blue() * v.blue();                        
        }    

        norms[c] = sqrt(magnitude);
    }
}

int LightTransportMatrix::get_cols() const {
    return m.get_cols(); 
}

int LightTransportMatrix::get_rows() const {
    return m.get_rows();
}

void LightTransportMatrix::resize(int rows, int cols) {
    m.resize(rows, cols);
}

void LightTransportMatrix::save(const char *file) const {
    m.save(file);
}

} // end namespace

