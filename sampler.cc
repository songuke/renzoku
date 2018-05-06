#include "sampler.h"
#include "vec2.h"
#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace std;

namespace Renzoku {

/**
 * SamplerIteratorImpl needs to access the internal of Sampler to retrieve samples stored in
 * a vector or list. In Java, implementations of Iterator interface are implemented as 
 * nested classes in the parent collection class (e.g., List). Here we declare 
 * SamplerIteratorImpl as a friend class of Sampler to achieve similar purpose.
 */
class SamplerIteratorImpl : public SamplerIterator {
protected:
    SamplerIteratorImpl(Sampler *s) : s(s) {}
    friend class Sampler;   // allow Sampler to access the above constructor

public:
    bool has_next_1D() {
        return s->sample_1D_index < s->sample_1Ds.size();
    }

    bool has_next_2D() {
        return s->sample_2D_index < s->sample_2Ds.size();        
    }

    Float next_1D() {
        if (has_next_1D())
            return s->sample_1Ds[s->sample_1D_index++];
        else
            throw std::out_of_range("Not enough 1D samples.");
    }

    Vec2 next_2D() {
        if (has_next_2D())
            return s->sample_2Ds[s->sample_2D_index++];
        else
            throw std::out_of_range("Not enough 2D samples.");
    }
    
protected:
    Sampler *s;
};

SamplerIterator *Sampler::get_iterator() {
    // FIXME: memory leak if forget to delete iterator after use. 
    return new SamplerIteratorImpl(this);
}

void Sampler::reset_1D() {
    sample_1Ds.clear();
    sample_1D_index = 0;
}

void Sampler::reset_2D() {
    sample_2Ds.clear();
    sample_2D_index = 0;
}

}

