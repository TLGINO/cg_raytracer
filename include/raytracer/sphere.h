#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "../../libs/glm/glm.hpp"
#include "object.h"

namespace raytracer {

/**
 * Implementation of the class Object for sphere shape.
 */
class Sphere : public Object {
 private:
  float radius_;
  vec3 center_;

 public:
  Sphere(float, vec3, vec3);
  Sphere(float, vec3, Material);
  Hit intersect(Ray) override;
};

}  // namespace raytracer

#endif  // RAYTRACER_SPHERE_H
