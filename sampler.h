#ifndef _SAMPLER_H_
#define _SAMPLER_H_

#include "common.h"
#include "random.h"

namespace Renzoku {

class SamplerIterator {    
public:
    virtual bool has_next_1D() = 0;
    virtual bool has_next_2D() = 0;
    virtual Float next_1D() = 0;
    virtual Vec2  next_2D() = 0;
};

/**
 * \brief Sampler is a container that stores a collection of samples requested by integrators. 
 * A sampler should be able to generate 1D and 2D samples. 
 * 
 * The integrator informs the sampler about the number of samples of each type that it needs at
 * the beginning of the rendering process. 
 *
 * In Mitsuba, Sampler class does not know what type of samples it generates except 1D or 2D samples. 
 * It is possible to let Sampler be a proxy that returns all types of samples to integrators, such as camera or light samples. 
 * This requires Sampler to have more prior knowledge about samples it generates, e.g., camera models (pinhole, thin lens), etc. 
 * I feel that this increases the coupling between Sampler and other classes. Therefore, I follow the approach taken by Mitsuba. 
 */
class Sampler {
public:
    Sampler(Random &rd) : rd(rd) {}
    
    virtual Sampler *clone() const = 0;

    /**
     * Add a set of new samples to the current sample list. Stratified sampling can be supported. 
     */
    virtual void generate_1D(int count) = 0; 
    virtual void generate_2D(int count) = 0;

    /**
     * Return a new iterator to the current sample list. 
     * Use after generate() is called to traverse through all the samples generated.
     */
    SamplerIterator *get_iterator();

    /** 
     * Generate a single sample independently.
     */
    virtual Float generate_1D() = 0;
    virtual Vec2  generate_2D() = 0;

protected:
    /**
     * Clear all existing samples and reset index to zero.
     */
    void reset_1D();
    void reset_2D();

protected:
    typedef vector<Float>        Sample1DList;
    typedef vector<Vec2>         Sample2DList;
    
    Sample1DList        sample_1Ds;
    Sample2DList        sample_2Ds;
    
    int sample_1D_index;
    int sample_2D_index;
    
    friend class SamplerIteratorImpl;
    Random &rd;
};

}; // end namespace
#endif

