#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H

#include "../../libs/glm/glm.hpp"
#include <cmath>

namespace raytracer::texture {

using glm::vec3, glm::vec2;

vec3 checkerboard(vec2);
vec3 rainbow(vec2);

}  // namespace raytracer::texture

#endif  //  RAYTRACER_TEXTURE_H
