//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h

#include "glm/glm.hpp"

using glm::vec3;
using glm::vec2;

vec3 checkerboardTexture(vec2 uv) {
    /* Exercise 2 (3 points) */
    float N = 40;
    float f = int(floor(N * uv.s) + floor(N * uv.t)) % 2;
    return vec3(f);
}

vec3 rainbowTexture(vec2 uv) {
    /* Exercise 2 (5 points) */
    float N = 40;
    int result = int(floor(N * uv.t + 0.5 * N * uv.s)) % 3;
    switch (result) {
    case 0:
        return {1.0f, 0.0f, 0.0f};
    case 1:
        return {0.0f, 1.0f, 0.0f};
    default:
        return {0.0f, 0.0f, 1.0f};
    }
}
#endif /* Textures_h */
