#ifndef RAYTRACER_LIGHTNING_H
#define RAYTRACER_LIGHTNING_H

#include "glm/glm.hpp"
#include "raytracer/light.h"
#include "raytracer/object.h"
#include "raytracer/ray.h"
#include <vector>

namespace raytracer::lightning {

using glm::vec3;
using std::vector;

vec3 toneMapping(vec3);
vec3 traceRay(vector<Object*>, Ray, vector<Light*>, vec3);

}  // namespace raytracer::lightning

#endif  // RAYTRACER_LIGHTNING_H
