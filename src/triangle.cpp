#include "../include/raytracer/triangle.h"

using glm::vec3;
using raytracer::Hit;
using raytracer::Triangle;

Triangle::Triangle(vec3 pa, vec3 pb, vec3 pc, Material mat)
    : point_a(pa), point_b(pb), point_c(pc) {
  this->setMaterial(mat);
}
Triangle::Triangle(vec3 pa, vec3 pb, vec3 pc, vec3 na, vec3 nb, vec3 nc,
                   Material mat)
    : point_a(pa),
      point_b(pb),
      point_c(pc),
      normal_a(na),
      normal_b(nb),
      normal_c(nc) {
  this->setMaterial(mat);
  this->with_normal = true;
}

Hit Triangle::intersect(Ray ray) {
  Hit hit;

  vec3 e1 = point_b - point_a;
  vec3 e2 = point_c - point_a;
  vec3 h = cross(ray.direction, e2);
  float a = glm::dot(e1, h);

  if (a > -EPSILON && a < EPSILON) return hit;

  float f = 1.0f / a;
  vec3 s = ray.origin - point_a;
  float u = f * dot(s, h);

  if (u < 0.0f || u > 1.0f) return hit;

  vec3 q = cross(s, e1);
  float v = f * dot(ray.direction, q);

  if (v < 0.0f || u + v > 1.0f) return hit;

  float t = f * dot(e2, q);
  if (t > EPSILON) {
    hit.hit = true;
    hit.intersection = ray.origin + t * ray.direction;
    if (this->with_normal) {
      vec3 cross_product = cross({point_b - point_a}, (point_c - point_a));
      float W = 0.5 * length(cross_product);

      vec3 center = (point_a + point_b + point_c) / vec3(3.0);
      float lambda_1 =
          (0.5 * (length(cross({point_b - center}, (point_c - center))))) / W;
      float lambda_2 =
          (0.5 * (length(cross({point_c - center}, (point_a - center))))) / W;
      // lambda3 received implicitly from difference of the other 2: 1 - l1 - l2

      hit.normal = normalize((1.0f - lambda_1 - lambda_2) * normal_a +
                             lambda_1 * normal_b + lambda_2 * normal_c);
    } else {
      hit.normal = normalize(cross(point_b - point_a, point_c - point_a));
    }
    hit.distance = distance(ray.origin, hit.intersection);
    hit.object = this;
    return hit;
  }
  return hit;
}
