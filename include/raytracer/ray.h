#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H

#include "../../libs/glm/glm.hpp"

namespace raytracer {

using glm::vec3;

/**
 * Class representing a single ray.
 */
class Ray {
 public:
  vec3 origin;
  vec3 direction;

  Ray(vec3, vec3);
};

}  // namespace raytracer

#endif  // RAYTRACER_RAY_H
