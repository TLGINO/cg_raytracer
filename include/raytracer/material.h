//
//  Material.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include "../../libs/glm/glm.hpp"

namespace raytracer {

using glm::vec2;
using glm::vec3;

/**
 * Structure describing a material of an object
 */
struct Material {
  vec3 ambient = vec3(0.0F);   // Ambient color of the object
  vec3 diffuse = vec3(1.0F);   // Color of the object
  vec3 specular = vec3(0.0F);  // Specular color of the object
  float shininess = 0.0F;      // Exponent for Phong model
  vec3 (*texture)(vec2 uv) = nullptr;
};

}  // namespace raytracer

#endif  // RAYTRACER_MATERIAL_H
