/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

#include "Image.h"
#include "Material.h"

using namespace std;
using glm::mat4;
using glm::vec3;
using glm::vec4;
using glm::normalize;

/**
 Class representing a single ray.
 */
class Ray
{
public:
  glm::vec3 origin;    ///< Origin of the ray
  glm::vec3 direction; ///< Direction of the ray
                       /**
                        Contructor of the ray
                        @param origin Origin of the ray
                        @param direction Direction of the ray
                        */
  Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction)
  {
  }
};

class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit
{
  bool hit;               ///< Boolean indicating whether there was or there was no intersection with an object
  glm::vec3 normal;       ///< Normal vector of the intersected object at the intersection point
  glm::vec3 intersection; ///< Point of Intersection
  float distance;         ///< Distance from the origin of the ray to the intersection point
  Object *object;         ///< A pointer to the intersected object
  glm::vec2 uv;           ///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object {
protected:
  glm::mat4 transformationMatrix;        ///< Matrix representing the transformation from the local to the global coordinate system
  glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
  glm::mat4 normalMatrix;                ///< Matrix for transforming normal vectors from the local to the global coordinate system

public:
  glm::vec3 color;   ///< Color of the object
  Material material; ///< Structure describing the material of the object
                     /** A function computing an intersection, which returns the structure Hit */
  virtual Hit intersect(Ray ray) = 0;

  /** Function that returns the material struct of the object*/
  Material getMaterial()
  {
    return material;
  }
  /** Function that set the material
   @param material A structure describing the material of the object
  */
  void setMaterial(Material material)
  {
    this->material = material;
  }
  /** Functions for setting up all the transformation matrices
   @param matrix The matrix representing the transformation of the object in the global coordinates */
  void setTransformation(mat4 matrix) {
    transformationMatrix = matrix;
    inverseTransformationMatrix = glm::inverse(transformationMatrix);
    normalMatrix = glm::transpose(inverseTransformationMatrix);
  }

  /**
   * Transform ray into local coords:
   * affine M-1 * homogenous coords of ray,
   * then cast back to vec3 and normalize direction
   * @param ray ray in global coords
   * @return ray in local coords
   */
  Ray toLocal(const Ray& ray) {
    return Ray(
        vec3(inverseTransformationMatrix * vec4(ray.origin, 1)),
        normalize(vec3(inverseTransformationMatrix * vec4(ray.direction, 0))));
  }

  /**
   * Transform hit into global coords:
   * M * intersection and M^(-t) * normal to give it back in global coords
   * @param hit hit in local coords
   * @return hit in global coords
   */
  Hit toGlobal(const Hit& h, const vec3& ray_origin) {
    return {
        .hit = h.hit,
        .normal = normalize(vec3(normalMatrix * vec4(h.normal, 0))),
        .intersection = vec3(transformationMatrix * vec4(h.intersection, 1)),
        .distance = glm::distance(h.intersection, ray_origin),
        .object = h.object,
        .uv = h.uv,
    };
  }
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object {
private:
  float radius = 1;                ///< Radius of the sphere
  vec3 center = vec3(0); ///< Center of the sphere

public:
  /**
   The constructor of the sphere
   @param color Color of the sphere
   @param material Material of the sphere
   */
  Sphere(vec3 color)
  {
    this->color = color;
  }
  Sphere(Material material)
  {
    this->material = material;
  }

  /** Implementation of the intersection function*/
  Hit intersect(Ray ray) {
    Ray local_ray = toLocal(ray);
    vec3 c = center - local_ray.origin;

    float cdotc = glm::dot(c, c);
    float cdotd = glm::dot(c, local_ray.direction);
    float D = 0;
    if (cdotc > cdotd * cdotd) {
      D = sqrt(cdotc - cdotd * cdotd);
    }

    if (D > radius) {
      return Hit();
    }

    float t1 = cdotd - sqrt(radius * radius - D * D);
    float t2 = cdotd + sqrt(radius * radius - D * D);
    float t = t1;
    if (t < 0)
      t = t2;

    if (t < 0) {
      return Hit();
    }

    Hit hit;
    hit.hit = true;
    hit.intersection = local_ray.origin + t * local_ray.direction;
    hit.normal = glm::normalize(hit.intersection - center);
    hit.distance = glm::distance(local_ray.origin, hit.intersection);
    hit.object = this;
    hit.uv.s = (asin(hit.normal.y) + M_PI / 2) / M_PI;
    hit.uv.t = (atan2(hit.normal.z, hit.normal.x) + M_PI) / (2 * M_PI);
    return toGlobal(hit, local_ray.origin);
  }
};

class Plane : public Object {
private:
  glm::vec3 normal;
  glm::vec3 point;

public:
  Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal) {}
  Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal)
  {
    this->material = material;
  }
  Hit intersect(Ray ray) {
    Hit hit;
    hit.hit = false;
    float DdotN = glm::dot(ray.direction, normal);
    if (DdotN < 0) {
      float PdotN = glm::dot(point - ray.origin, normal);
      float t = PdotN / DdotN;

      if (t > 0) {
        hit.hit = true;
        hit.normal = normal;
        hit.distance = t;
        hit.object = this;
        hit.intersection = t * ray.direction + ray.origin;
      }
    }
    return hit;
  }
};

class Cone : public Object
{
private:
  Plane *cover;

public:
  Cone(Material material)
  {
    this->material = material;
    this->cover = new Plane(glm::vec3(0, 1, 0), glm::vec3(0, 1, 0));
  }
  Hit intersect(Ray ray)
  {

    Hit hit;
    hit.hit = false;

    /*  ---- Assignemnt 5 -----

     Implement the ray-cone intersection. Before intersecting the ray with the cone,
     make sure that you transform the ray into the local coordinate system.
     Remember about normalizing all the directions after transformations.

    */

    /* If the intersection is found, you have to set all the critical fields in the Hit strucutre
     Remember that the final information about intersection point, normal vector and distance have to be given
     in the global coordinate system.
     */

    glm::vec3 n_direction = this->inverseTransformationMatrix * glm::vec4(ray.direction, 0.0);
    glm::vec3 n_origin = this->inverseTransformationMatrix * glm::vec4(ray.origin, 1.0);
    n_direction = normalize(n_direction);

    float a = (n_direction.x * n_direction.x) +
              (n_direction.z * n_direction.z) -
              (n_direction.y * n_direction.y);

    float b = 2 *
              ((n_direction.x * n_origin.x) + (n_direction.z * n_origin.z) - (n_direction.y * n_origin.y));

    float c = (n_origin.x * n_origin.x) +
              (n_origin.z * n_origin.z) -
              (n_origin.y * n_origin.y);

    float delta = pow(b, 2) - (4 * a * c);

    if (delta < 0)
      return hit;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    hit.intersection = n_origin + t1 * n_direction;
    if (t1 < 0 || hit.intersection.y > 1 || hit.intersection.y < 0)
    {
      hit.intersection = n_origin + t2 * n_direction;
      if (hit.intersection.y > 1 || hit.intersection.y < 0 || t2 < 0)
        return hit;
    };

    hit.normal = glm::vec3(hit.intersection.x, -hit.intersection.y, hit.intersection.z);
    hit.normal = normalize(hit.normal);

    Ray new_ray(n_origin, n_direction);
    Hit hit_plane = this->cover->intersect(new_ray);
    if (hit_plane.hit && hit_plane.distance < t1 && length(hit_plane.intersection - glm::vec3(0, 1, 0)) <= 1.0)
    {
      hit.intersection = hit_plane.intersection;
      hit.normal = hit_plane.normal;
    }

    hit.hit = true;
    hit.object = this;
    hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1.0);
    hit.normal = (normalMatrix * glm::vec4(hit.normal, 0.0));
    hit.normal = normalize(hit.normal);
    hit.distance = length(hit.intersection - ray.origin);

    return hit;
  }
};

/**
 Light class
 */
class Light
{
public:
  glm::vec3 position; ///< Position of the light source
  glm::vec3 color;    ///< Color/intentisty of the light source
  Light(glm::vec3 position) : position(position)
  {
    color = glm::vec3(1.0);
  }
  Light(glm::vec3 position, glm::vec3 color) : position(position), color(color)
  {
  }
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001, 0.001, 0.001);
vector<Object *> objects; ///< A list of all objects in the scene

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction, Material material)
{

  glm::vec3 color(0.0);
  for (int light_num = 0; light_num < lights.size(); light_num++)
  {

    glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
    glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

    float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
    float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

    glm::vec3 diffuse_color = material.diffuse;
    if (material.texture)
    {
      diffuse_color = material.texture(uv);
    }

    glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
    glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));

    // distance to the light
    float r = glm::distance(point, lights[light_num]->position);
    r = max(r, 0.1f);

    color += lights[light_num]->color * (diffuse + specular) / r / r;
  }
  color += ambient_light * material.ambient;

  return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray)
{

  Hit closest_hit;

  closest_hit.hit = false;
  closest_hit.distance = INFINITY;

  for (int k = 0; k < objects.size(); k++)
  {
    Hit hit = objects[k]->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance)
      closest_hit = hit;
  }

  glm::vec3 color(0.0);

  if (closest_hit.hit)
  {
    color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
  }
  else
  {
    color = glm::vec3(0.0, 0.0, 0.0);
  }
  return color;
}
/**
 Function defining the scene
 */
void sceneDefinition(float up_down_val, float rotation_val) {
  Material green_diffuse;//
  green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);//
  green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);//

  Material red_specular;//
  red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);//
  red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);//
  red_specular.specular = glm::vec3(0.5);//
  red_specular.shininess = 10.0;//

  Material blue_specular;//
  blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);//
  blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);//
  blue_specular.specular = glm::vec3(0.6);//
  blue_specular.shininess = 100.0;//

  Sphere *s1 = new Sphere(blue_specular);
  s1->setTransformation(translate(glm::vec3(1, -2, 8)));
  objects.push_back(s1);

  Sphere *s2 = new Sphere(red_specular);
  s2->setTransformation(
          translate(vec3(-0.7, -1, 3.5)) *
                scale(mat4(1), vec3(0.5))
                );
//      scale(
//          translate(mat4(1.0f), vec3(-0.20, -0.20, 0.6)),
//          glm::vec3(0.1)));
  objects.push_back(s2);


  // ------ Assignment 5 -------
  // You can remove the green sphere as it should be replaced with a green cone
  // objects.push_back(new Sphere(1.0, glm::vec3(3, -2, 6), green_diffuse));

  // Textured sphere
  Material textured;
  textured.texture = &rainbowTexture;
  Sphere *s3 = new Sphere(textured);
  s3->setTransformation(
//      translate(vec3(-2, 0, 35)) *
//            scale(mat4(1), vec3(10))
//            );
      // translate(glm::vec3(-6, 4, 23)) *


      translate(glm::vec3(-6, up_down_val, 23)) *
      scale(glm::vec3(7)) *
      rotate(glm::radians(rotation_val), glm::vec3(1, 0, 0)));
  // rotate(glm::radians(90.0f), glm::vec3(1, 0, 0)));


  objects.push_back(s3);

  // Planes
  Material red_diffuse;//
  red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);//
  red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);//

  Material blue_diffuse;//
  blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);//
  blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);//
  objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0.0, 1, 0)));//
  objects.push_back(new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), green_diffuse));//
  objects.push_back(new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));//
  objects.push_back(new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));//
  objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));//
  objects.push_back(new Plane(glm::vec3(0, 1, -0.01), glm::vec3(0.0, 0.0, 1.0), green_diffuse));//

  /* ----- Assignment 5 -------
  Create two cons and add them to the collection of our objects.
  Remember to create them with corresponding materials and transformation matrices
  */
  Material yellow_very_specular;
  yellow_very_specular.ambient = glm::vec3(0.1f, 0.10f, 0.0f);
  yellow_very_specular.diffuse = glm::vec3(0.5f, 0.5f, 0.0f);
  yellow_very_specular.specular = glm::vec3(1.0);
  yellow_very_specular.shininess = 100.0;

  Cone *cone1 = new Cone(yellow_very_specular);
  cone1->setTransformation(translate(glm::vec3(5, 9, 14)) *
                           scale(glm::vec3(3.0f, 12.0f, 3.0f)) *
                           rotate(glm::radians(180.0f), glm::vec3(1, 0, 0)));
  objects.push_back(cone1);

  Cone *cone2 = new Cone(green_diffuse);
  cone2->setTransformation(translate(glm::vec3(6, -3, 7)) *
                           rotate(glm::atan(3.0f), glm::vec3(0, 0, 1)) *
                           scale(glm::vec3(1.0f, 3.0f, 1.0f)));
  objects.push_back(cone2);

  lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));//
  lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));//
  lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));//
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity)
{
  float gamma = 1.0 / 2.0;
  float alpha = 12.0f;
  return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char *argv[])
{

  clock_t t = clock(); // variable for keeping the time of the rendering

  int width = 1024; // width of the image
  int height = 768; // height of the image
  float fov = 90;   // field of view
  if (argc == 4)
  {
    sceneDefinition(stoi(argv[2]), stoi(argv[3])); // Let's define a scene
  }
  else
  {

    sceneDefinition(4, 0); // Let's define a scene
  }

  Image image(width, height); // Create an image where we will store the result

  float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
  float X = -s * width / 2;
  float Y = s * height / 2;

  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++)
    {

      float dx = X + i * s + s / 2;
      float dy = Y - j * s - s / 2;
      float dz = 1;

      glm::vec3 origin(0, 0, 0);
      glm::vec3 direction(dx, dy, dz);
      direction = glm::normalize(direction);

      Ray ray(origin, direction);

      image.setPixel(i, j, toneMapping(trace_ray(ray)));
    }

  t = clock() - t;
  cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
  cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

  // Writing the final results of the rendering
  if (argc >= 2)
  {
    string res = argv[1];
    res = res + ".ppm";
    image.writeImage(res.c_str());
  }
  else
  {
    image.writeImage("./result.ppm");
  }

  return 0;
}
