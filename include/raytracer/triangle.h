#ifndef RAYTRACER_TRIANGLE_H
#define RAYTRACER_TRIANGLE_H

#include "../../libs/glm/glm.hpp"
#include "hit.h"
#include "material.h"
#include "object.h"

namespace raytracer {

class Triangle : public Object {
 public:
  Triangle(vec3, vec3, vec3, Material);
  Triangle(vec3, vec3, vec3, vec3, vec3, vec3, Material);
  Hit intersect(Ray) override;

 private:
  vec3 point_a, point_b, point_c;
  vec3 normal_a = vec3(0), normal_b = vec3(0), normal_c = vec3(0);
  float const EPSILON = 1e-6;
  bool with_normal = false;
};

}  // namespace raytracer

#endif  // RAYTRACER_TRIANGLE_H
