#include "image_sampler.h"
#include "scene.h"
#include "camera.h"
#include "random.h"

namespace Renzoku {

void ImageSampler::sample_stratified(Scene *scene, int block_size, Pixels &pixels) {
    Camera *camera = scene->get_camera();
    Random &rd     = *scene->get_random();

    Size2 size = camera->get_image_size();
    
    for (int i = 0; i < size.height; i += block_size) {
        for (int j = 0; j < size.width; j += block_size) {
            Vec2 sample;
            sample.random(rd);
            int x = j + sample.x() * block_size;
            int y = i + sample.y() * block_size;

            pixels.push_back(Vec2(x, y));            
        }
    }
}

} // end namespace
