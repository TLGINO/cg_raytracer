#include "../include/raytracer/light.h"

using raytracer::Light;

Light::Light(vec3 position, vec3 color) : position(position), color(color) {}
Light::Light(vec3 position) : Light(position, vec3(1.0F)) {}
