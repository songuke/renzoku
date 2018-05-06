#include "light_particle.h"
#include "light_tree.h"

namespace Renzoku {

LightParticle::LightParticle() {
    p = Vec3();
}

LightParticle::LightParticle(Photon *photon) : photon(photon) {
    p = photon->p;
}

LightParticle::LightParticle(BrdfPointLight *light) : light(light) {
    p = light->org();
}

LightParticle::LightParticle(DensityPoint *point) : point(point) {
    p = point->p;
}

LightParticle::LightParticle(PathNode *node) : node(node) {
    p = node->dg.p;
}

LightParticle::LightParticle(LightNode *node) : light_node(node) {
    p = light_node->bb.centroid();
}

} // end namespace Renzoku
