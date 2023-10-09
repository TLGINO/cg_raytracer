//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h

#include "glm/glm.hpp"

glm::vec3 rainbowTexture(glm::vec2 uv)
{
    float n = 20;
    int value = int(floor(n * uv.t + 1.5 * n * uv.s)) % 3;
    switch (value)
    {
    case 0:
        return glm::vec3(0, 0, 1);
    case 1:
        return glm::vec3(0, 1, 0);
    default:
        return glm::vec3(1, 0, 0);
    }
}
#endif /* Textures_h */
