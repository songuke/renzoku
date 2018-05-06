#ifndef _LIGHT_TRANSPORT_H_
#define _LIGHT_TRANSPORT_H_

#include "common.h"

#include "matrix.h"
#include "matrix_c.h"

#include "vpl_interface.h"

#include "image_sampler.h"

namespace Renzoku {

/**
 * A common interface for light transport matrix and its components like form factor and BRDF matrix.
 */
class TransportMatrix {
public:
    /**
     * Generate light transport matrix content from a set of pixel samples and a set of VPLs.
     */
    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels) = 0;
    
    /**
     * Return the distance between two columns
     */
    virtual Float get_distance(int i, int j) { return 0.; }
    
    /**
     * Generate the probability distribution for clustering by sampling. 
     *
     * Alpha vector is calculated by Hasan's SG07 approach.
     *
     * This function is written here for better performance. 
     */
    virtual void generate_alphas_for_clustering(vector<Float> &alphas) {}

    /**
     * Return the L2-norm of the column
     */
    virtual Float get_magnitude_of_column(int i) { return 0.; }

    virtual int get_cols() const = 0;
    virtual int get_rows() const = 0;

    virtual void save(const char *file) const = 0;
};

/**
 * Light transport is a RGB matrix that stores the out-going radiance 
 * of each surface point due to each virtual point light. 
 */
class LightTransportMatrix : public TransportMatrix {
public:    
    LightTransportMatrix(int rows, int cols);

    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels);
    
    virtual void resize(int rows, int cols);
    
    virtual void get_column(int i, vector<Float> &c) const;
    virtual void set_column(int i, const vector<Float> &c);

    virtual void generate_alphas_for_clustering(vector<Float> &alphas);

    virtual Float get_column_norm(int i) const;
    virtual void get_column_norms(vector<Float> &norms) const;

    virtual int get_rows() const;
    virtual int get_cols() const;     
    virtual MatrixC<Rgb> &get_matrix() { return m; }

    virtual void save(const char *file) const;

protected:
    /**
     * Due to the frequent use of calculating distance between two columns,
     * a column-major internal storage is prefered for this matrix.
     */
    MatrixC<Rgb> m;
};


/**
 * Form factor matrix stores form factor values between each virtual point light and each surface.
 * No visibility tracing is done. 
 */ 
class FormFactorMatrix : public TransportMatrix {
public:    
    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels);
    
    virtual void generate_alphas_for_clustering(vector<Float> &alphas);

    virtual int get_rows() const;
    virtual int get_cols() const; 

    virtual void save(const char *file) const;

protected:
    MatrixC<Float> m;
};


/**
 * Visibility matrix is a 0/1 matrix that denotes the visiblity between each virtual point light 
 * and each surface point.
 *
 * TODO:
 * Derived from FormFactorMatrix class for code re-use. In fact, to save more storage, 
 * it can be derived directly from TransportMatrix and use a MatrixC<char> storage. 
 */ 
class VisibilityMatrix : public FormFactorMatrix {
public:    
    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels);
};


/**
 * This matrix captures all BRDF values at every virtual point light 
 * that contributes to every surface.
 */
class BrdfLightMatrix : public LightTransportMatrix {
public:    
    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels);    
};


/*
 * This matrix captures all BRDF values at every surface that receives 
 * contribution from every virtual point light. 
 */
class BrdfSurfaceMatrix : public LightTransportMatrix {
public:    
    virtual void fill(Scene *scene, IVirtualPointLightEvaluator *renderer, 
                      BrdfPointLights &vpls, Pixels &pixels);    
};


/**
 * A unified light transport class. 
 */
class LightTransport {
protected:
    LightTransportMatrix    T;
    FormFactorMatrix        G;
    BrdfLightMatrix         F1;
    BrdfSurfaceMatrix       F2;
};

} // end namespace

#endif