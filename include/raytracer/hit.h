#ifndef RAYTRACER_HIT_H
#define RAYTRACER_HIT_H

#include "../../libs/glm/glm.hpp"

namespace raytracer {

using glm::vec2;
using glm::vec3;

class Object;

/**
 * Structure representing the event of hitting an object
 */
struct Hit {
  bool hit = false;  // state if intersection with an object
  vec3 normal;  // Normal vector of intersected object at the intersection point
  vec3 intersection;          // Point of Intersection
  float distance = INFINITY;  // Distance from ray origin to intersection point
  Object* object = nullptr;   // pointer to the intersected object
  vec2 uv;  // Coordinates for computing the texture (texture coordinates)
};

}  // namespace raytracer

#endif  // RAYTRACER_HIT_H
