#include "../include/raytracer/ray.h"

using raytracer::Ray;

/**
 * Constructor of the ray
 * @param origin Origin of the ray
 * @param direction Direction of the ray
 */
Ray::Ray(vec3 origin, vec3 direction) : origin(origin), direction(direction) {}
