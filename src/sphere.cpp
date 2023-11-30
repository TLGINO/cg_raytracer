#include "../include/raytracer/sphere.h"

using raytracer::Hit;
using raytracer::Sphere;

/**
 * Constructor of the sphere
 * @param radius Radius of the sphere
 * @param center Center of the sphere
 * @param color Color of the sphere
 */
Sphere::Sphere(float radius, vec3 center, vec3 color)
    : radius_(radius), center_(center) {
  this->color = color;
}

Sphere::Sphere(float radius, vec3 center, Material material)
    : radius_(radius), center_(center) {
  this->material = material;
}

/**
 * Implementation of the intersection function
 */
Hit Sphere::intersect(Ray ray) {
  vec3 c = center_ - ray.origin;
  float cdotc = glm::dot(c, c);
  float cdotd = glm::dot(c, ray.direction);
  Hit hit;
  float D = 0;
  if (cdotc > cdotd * cdotd) {
    D = glm::sqrt(cdotc - cdotd * cdotd);
  }

  if (D > radius_) {
    hit.hit = false;
    return hit;
  }

  hit.hit = true;
  float t1 = cdotd - sqrt(radius_ * radius_ - D * D);
  float t2 = cdotd + sqrt(radius_ * radius_ - D * D);

  float t = t1;
  if (t < 0) t = t2;

  if (t < 0) {
    hit.hit = false;
    return hit;
  }

  hit.intersection = ray.origin + t * ray.direction;
  hit.normal = normalize(hit.intersection - center_);
  hit.distance = distance(ray.origin, hit.intersection);
  hit.object = this;

  // Ex2: computing texture coordinates for the sphere.
  hit.uv.s = (asin(hit.normal.y) + M_PI / 2) / M_PI;
  hit.uv.t = (atan2(hit.normal.z, hit.normal.x) + M_PI) / (2 * M_PI);
  return hit;
}
