/**
 * @file Image.h
 */

//  Image.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Image_h
#define Image_h

#include "glm/glm.hpp"
#include <fstream>

/**
 * Class allowing for creating an image and writing it to a file
 */
class Image {
 private:
  int width_, height_;  // width and height of the image
  int* data_;           // a pointer to the data representing the images

 public:
  /**
   * @param width with of the image
   * @param height height of the image
   */
  Image(int width, int height) : width_(width), height_(height) {
    data_ = new int[3 * width * height];
  }

  /**
   * Writes and image to a file in ppm format.
   * @param path the path where to the target image
   */
  void writeImage(char const* path) {
    std::ofstream file;
    file.open(path);
    file << "P3 \n";
    file << width_ << " " << height_ << '\n';
    file << 255 << '\n';
    for (int h = 0; h < height_; h++) {
      for (int w = 0; w < width_; w++) {
        file << data_[3 * (h * width_ + w)] << " ";
        file << data_[3 * (h * width_ + w) + 1] << " ";
        file << data_[3 * (h * width_ + w) + 2] << "  ";
      }
      file << '\n';
    }
    file.close();
  }

  /**
   * Set a value for one pixel
   * @param x x coordinate of the pixel - index of the column counting from left
   * to right
   * @param y y coordinate of the pixel - index of the row counting from top to
   * bottom
   * @param r red chanel value in range from 0 to 255
   * @param g green chanel value in range from 0 to 255
   * @param b blue chanel value in range from 0 to 255
   */
  [[maybe_unused]] void setPixel(int x, int y, int r, int g, int b) {
    data_[3 * (y * width_ + x)] = r;
    data_[3 * (y * width_ + x) + 1] = g;
    data_[3 * (y * width_ + x) + 2] = b;
  }

  /**
   * Set a value for one pixel
   * @param x x coordinate of the pixel - index of the column counting from left
   * to right
   * @param y y coordinate of the pixel - index of the row counting from top to
   * bottom
   * @param r red chanel value in range from 0 to 1
   * @param g green chanel value in range from 0 to 1
   * @param b blue chanel value in range from 0 to 1
   */
  [[maybe_unused]] void setPixel(int x, int y, float r, float g, float b) {
    data_[3 * (y * width_ + x)] = 255 * r;
    data_[3 * (y * width_ + x) + 1] = 255 * g;
    data_[3 * (y * width_ + x) + 2] = 255 * b;
  }

  /**
   * Set a value for one pixel
   * @param x x coordinate of the pixel - index of the column counting from left
   * to right
   * @param y y coordinate of the pixel - index of the row counting from top to
   * bottom
   * @param color color of the pixel expressed as vec3 of RGB values in range
   * from 0 to 1
   */
  void setPixel(int x, int y, glm::vec3 color) {
    data_[3 * (y * width_ + x)] = 255 * color.r;
    data_[3 * (y * width_ + x) + 1] = 255 * color.g;
    data_[3 * (y * width_ + x) + 2] = 255 * color.b;
  }
};

#endif /* Image_h */
