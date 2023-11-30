#ifndef RAYTRACER_OBJECT_H
#define RAYTRACER_OBJECT_H

#include "../../libs/glm/glm.hpp"
#include "hit.h"
#include "material.h"
#include "ray.h"

namespace raytracer {

/**
 * General class for the object
 */
class Object {
 public:
  vec3 color = vec3(0.0F);
  Material material;  // Structure describing the material of the object

  /**
   * pure virtual function computing intersection, returns the structure Hit
   */
  virtual Hit intersect(Ray) = 0;

  Material getMaterial() const;
  void setMaterial(Material);
};

}  // namespace raytracer

#endif  // RAYTRACER_OBJECT_H
