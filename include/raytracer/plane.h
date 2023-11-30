#ifndef RAYTRACER_PLANE_H
#define RAYTRACER_PLANE_H

#include "../../libs/glm/glm.hpp"
#include "hit.h"
#include "material.h"
#include "object.h"

namespace raytracer {

class Plane : public Object {
 private:
  vec3 point_ = vec3(0);
  vec3 normal_ = vec3(0, 0, 1);

 public:
  Plane(vec3, vec3);
  Plane(vec3, vec3, Material);
  Hit intersect(Ray) override;
};

}  // namespace raytracer

#endif  // RAYTRACER_PLANE_H
