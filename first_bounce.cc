#include "first_bounce.h"

namespace Renzoku {

FirstBounce::FirstBounce() { 
    suffix = "1st";
}

void FirstBounce::initialize(Scene *scene) {
    MonteCarloIntegrator::initialize(scene);
    //sampler = new StratifiedSampler(*scene->get_random());
}

}
