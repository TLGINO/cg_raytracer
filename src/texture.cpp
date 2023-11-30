#include "../include/raytracer/texture.h"

[[maybe_unused]] glm::vec3 checkerboard(glm::vec2 uv) {
  /* Exercise 2 (3 points) */
  float n = 40;
  float x = std::floor(n * uv.s);
  float y = std::floor(n * uv.t);
  float f = static_cast<int>(x + y) % 2;
  return glm::vec3(f);
}

glm::vec3 rainbow(glm::vec2 uv) {
  /* Exercise 2 (5 points) */
  float n = 40;
  int result = static_cast<int>(std::floor(n * uv.t + 0.5 * n * uv.s)) % 3;
  switch (result) {
    case 0:
      return {1.0F, 0.0F, 0.0F};
    case 1:
      return {0.0F, 1.0F, 0.0F};
    default:
      return {0.0F, 0.0F, 1.0F};
  }
}
