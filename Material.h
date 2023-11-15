//
//  Material.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Material_h
#define Material_h

#include "glm/glm.hpp"
#include "Textures.h"

/**
 * Structure describing a material of an object
 */
struct Material
{
    vec3 ambient = vec3(0.0f);  // Ambient color of the object
    vec3 diffuse = vec3(1.0f);  // Color of the object
    vec3 specular = vec3(0.0f); // Specular color of the object
    float shininess = 0.0f;     // Exponent for Phong model
    vec3 (*texture)(vec2 uv) = nullptr;
    float is_reflective = 0.0f;
    float refraction_index = false;
    float refraction = false;
};

#endif /* Material_h */
