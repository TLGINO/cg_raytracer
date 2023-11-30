#ifndef RAYTRACER_IMAGE_H
#define RAYTRACER_IMAGE_H

#include "glm/glm.hpp"
#include <fstream>

namespace raytracer {

/**
 * Class allowing for creating an image and writing it to a file
 */
class Image {
 private:
  [[maybe_unused]] int width_, height_;  // width and height of the image
  [[maybe_unused]] int* data_;  // a pointer to the data representing the images

 public:
  Image(int, int);

  void writeImage(char const*);
  [[maybe_unused]] void setPixel(int, int, int, int, int);
  [[maybe_unused]] void setPixel(int, int, float, float, float);
  void setPixel(int, int, glm::vec3);
};

}  // namespace raytracer

#endif  // RAYTRACER_IMAGE_H
