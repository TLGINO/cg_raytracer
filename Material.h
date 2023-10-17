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
struct Material {
    vec3 ambient = vec3(0.0);
    vec3 diffuse = vec3(1.0);
    vec3 specular = vec3(0.0);
    float shininess = 0.0; // Exponent for Phong model
    vec3 (* texture)(vec2 uv) = nullptr;
};

#endif /* Material_h */
