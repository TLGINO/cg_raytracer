#ifndef RAYTRACER_LIGHT_H
#define RAYTRACER_LIGHT_H

#include "../../libs/glm/glm.hpp"

namespace raytracer {

using glm::vec3;

/**
 * Light class
 */
class Light {
 public:
  vec3 position;
  vec3 color;

  Light(vec3, vec3);
  explicit Light(vec3);
};

}  // namespace raytracer

#endif  // RAYTRACER_LIGHT_H
