#include "sampler.h"
#include "vec2.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

namespace Renzoku {

class StratifiedSampler : public Sampler {
public:
    StratifiedSampler(Random &rd) : Sampler(rd) {}

    Sampler *clone() const {
        return new StratifiedSampler(rd);
    }

    void generate_1D(int count) {
        reset_1D();
        
        // TODO: permute the order of strata so we don't always get the sample from the first strata first.
        // for each strata, generate a sample within the strata
        for (int i = 0; i < count; ++i) {             
            Float v = (i + rd()) / count;
            sample_1Ds.push_back(v);
        }
    }

    void generate_2D(int count) {
        reset_2D();
                
        int size = ceil(sqrt(count));
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Float x = (j + rd()) / size;
                Float y = (i + rd()) / size;
                sample_2Ds.push_back(Vec2(x, y));
            }
        }
    }

    Float generate_1D() {
        throw std::runtime_error("Not implemented for stratified sampling.");
    }

    Vec2 generate_2D() {
        throw std::runtime_error("Not implemented for stratified sampling.");
    }
};

}
