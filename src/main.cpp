#include "glm/glm.hpp"
#include "raytracer/image.h"
#include "raytracer/light.h"
#include "raytracer/lightning.h"
#include "raytracer/mesh.h"
#include "raytracer/texture.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

#define PI atan(1) * 4

using glm::normalize;
using glm::vec3;
using raytracer::Image;
using raytracer::Light;
using raytracer::Material;
using raytracer::Mesh;
using raytracer::Object;
using raytracer::Ray;
using raytracer::lightning::toneMapping;
using raytracer::lightning::traceRay;
using raytracer::texture::rainbow;
using std::vector;

// GLOBAL VARIABLES
vector<Light*> lights;  // list of lights in the scene
vec3 ambient_light(0.5F);
vector<Object*> objects;  // list of all objects in the scene

/**
 * Function defining the scene
 */
void sceneDefinition() {
  vec3 color_red{1.5f, 0.3f, 0.3f};
  vec3 color_blue{0.4f, 0.4f, 1.5f};
  vec3 color_green{0.5f, 1.5f, 0.5f};
  vec3 color_white{1.5f, 1.5f, 1.5f};

  //  Material red_specular{
  //      .ambient = {0.01f, 0.03f, 0.03f},
  //      .diffuse = color_red,
  //      .specular = vec3(0.5f),
  //      .shininess = 10.0f,
  //  };

  Material blue_shiny{
      .ambient = {0.07f, 0.07f, 0.1f},
      .diffuse = color_blue,
      .specular = vec3(0.6f),
      .shininess = 100.0f,
  };

  //  Material green_diffuse{
  //      .ambient = {0.07f, 0.09f, 0.07f},
  //      .diffuse = color_green,
  //  };
  //
  //  Material white_plain{
  //      .ambient = {1.0f, 1.0f, 1.0f},
  //      .diffuse = color_white,
  //      .specular = vec3(1.6f),
  //      .shininess = 0.0f,
  //  };
  //  Material material_rainbow{
  //      .texture = &rainbow,
  //  };

  // TEST
  objects.push_back(new Mesh("./meshes/bunny_with_normals.obj", blue_shiny));
  // objects.push_back(new Mesh("./meshes/bunny.obj", blue_shiny));

  // objects.push_back(new Mesh("./meshes/armadillo.obj", blue_shiny));
  // objects.push_back(new Mesh("./meshes/lucy.obj", green_diffuse));
  // objects.push_back(new Triangle(
  // 	{-1.0f, -2.5f, 6.0f},
  // 	{1.0f, -2.0f, 8.0f},
  // 	{3.0f, -2.0f, 6.0f},
  // 	blue_shiny));

  // objects.push_back(new Triangle(
  // 	{-1.0f, 2.5f, 6.0f},
  // 	{1.0f, 2.0f, 8.0f},
  // 	{3.0f, 2.0f, 6.0f},
  // 	red_specular));

  // objects.push_back(new Triangle(
  // 	{-1.0f, 5.0f, 6.0f},
  // 	{-10.0f, 25.0f, 10.0f},
  // 	{-3.0f, 1.0f, 6.0f},
  // 	green_diffuse));

  // spheres
  // objects.push_back(new Sphere(0.5f, {-1.0f, -2.5f, 6.0f}, white_plain));
  // objects.push_back(new Sphere(0.5f, {-1.0f, -2.5f, 6.0f}, red_specular));
  // objects.push_back(new Sphere(1.0f, {1.0f, -2.0f, 8.0f}, blue_shiny));
  // objects.push_back(new Sphere(1.0f, {3.0f, -2.0f, 6.0f}, green_diffuse));
  // objects.push_back(new Sphere(6.0f, {-5.0f, 3.5f, 20.0f},
  // material_rainbow));

  // planes
  // objects.push_back(new Plane({0.0f, 0.0f, -0.01f},
  // 							{0.0f, 0.0f, 1.0f},
  // 							Material())); // back
  // objects.push_back(new Plane({15.0f, 0.0f, 0.0f},
  // 							{-1.0f, 0.0f, 0.0f},
  // 							{
  // 								.diffuse =
  // {0.6f, 0.6f, 1.0f},
  // 							})); // right
  // objects.push_back(new Plane({0.0f, 0.0f, 30.0f},
  // 							{0.0f, 0.0f, -1.0f},
  // 							{
  // 								.diffuse =
  // {0.5f, 1.0f, 0.5f},
  // 							})); // front
  // objects.push_back(new Plane({-15.0f, 0.0f, 0.0f},
  // 							{1.0f, 0.0f, 0.0f},
  // 							{
  // 								.diffuse =
  // {0.6f, 0.4f, 0.4f},
  // 							})); // left
  // objects.push_back(new Plane({0.0f, 27.0f, 0.0f},
  // 							{0.0f, -1.0f, 0.0f},
  // 							Material())); // top
  // objects.push_back(new Plane({0.0f, -3.0f, 0.0f},
  // 							{0.0f, 1.0f, 0.0f},
  // 							Material())); // bottom

  // lights
  lights.push_back(new Light({0.0F, 26.0F, 5.0F}, vec3(150.0F)));
  lights.push_back(new Light({0.0F, 1.0F, 12.0F}, vec3(25.0F)));
  lights.push_back(new Light({0.0F, 5.0F, 1.0F}, vec3(40.0F)));
}

int main(int argc, char const* argv[]) {
  clock_t t = clock();  // keeping the time of the rendering
  int const width = 1024;
  int const height = 768;
  float const fov = 90;  // field of view
  // Debug
  /*
    int width = 512;
    int height = 384;  // Final
  */
  // int width = 2048;
  // int height = 1536;

  sceneDefinition();
  Image image(width, height);  // Create an image where we will store the result
  float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
  float X = -s * width / 2;
  float Y = s * height / 2;

  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++) {
      float dx = X + i * s + s / 2;
      float dy = Y - j * s - s / 2;
      // vec3 origin(0);
      vec3 origin(0, 1, -4);
      vec3 direction(dx, dy, 1);
      direction = normalize(direction);
      Ray ray(origin, direction);
      image.setPixel(
          i, j, toneMapping(traceRay(objects, ray, lights, ambient_light)));
    }

  t = static_cast<float>(clock() - t);
  std::cout << "It took " << t / CLOCKS_PER_SEC
            << " seconds to render the image.\n";
  std::cout << "I could render at " << CLOCKS_PER_SEC / t
            << " frames per second.\n";

  // Writing the final results of the rendering
  if (argc == 2)
    image.writeImage(argv[1]);
  else
    image.writeImage("./result.ppm");

  return 0;
}
