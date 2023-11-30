#include "../include/raytracer/plane.h"

using raytracer::Hit;
using raytracer::Plane;

Plane::Plane(vec3 point, vec3 normal) : point_(point), normal_(normal) {}

Plane::Plane(vec3 point, vec3 normal, Material material)
    : point_(point), normal_(normal) {
  this->material = material;
}

/**
 * Ex1: Plane-ray intersection
 */
Hit Plane::intersect(Ray ray) {
  float ray_dot_n = dot(point_ - ray.origin, normal_);
  float d_dot_N = dot(ray.direction, normal_);

  // cos b of angle parallel to plane
  if (d_dot_N == 0) return Hit();

  float t = ray_dot_n / d_dot_N;

  // intersection behind ray origin
  if (t < 0) return Hit();

  Hit hit{
      .hit = true,
      .normal = normalize(normal_),
      .intersection = ray.origin + t * ray.direction,
      .distance = distance(ray.origin, hit.intersection),
      .object = this,
      .uv = {0, 0},
  };
  return hit;
}
