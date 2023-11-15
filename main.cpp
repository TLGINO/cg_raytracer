/**
@file main.cpp
*/

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "glm/gtx/vector_angle.hpp"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "Image.h"
#include "Material.h"

using namespace std;

/**
 Class representing a single ray.
 */
class Ray
{
public:
  glm::vec3 origin;    ///< Origin of the ray
  glm::vec3 direction; ///< Direction of the ray unit vector
                       /**
                        Contructor of the ray
                        @param origin Origin of the ray
                        @param direction Direction of the ray
                        */
  Ray(glm::vec3 origin, glm::vec3 direction)
      : origin(origin), direction(direction) {}
};

class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit
{
  bool hit;               ///< Boolean indicating whether there was or there was no
                          ///< intersection with an object
  glm::vec3 normal;       ///< Normal vector of the intersected object at the
                          ///< intersection point
  glm::vec3 intersection; ///< Point of Intersection
  float distance;         ///< Distance from the origin of the ray to the intersection
                          ///< point
  Object *object;         ///< A pointer to the intersected object
  glm::vec2 uv;           ///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object
{

protected:
  glm::mat4
      transformationMatrix;              ///< Matrix representing the transformation from the
                                         ///< local to the global coordinate system
  glm::mat4 inverseTransformationMatrix; ///< Matrix representing the
                                         ///< transformation from the global to
                                         ///< the local coordinate system
  glm::mat4 normalMatrix;                ///< Matrix for transforming normal vectors from the
                                         ///< local to the global coordinate system

public:
  glm::vec3 color;   ///< Color of the object
  Material material; ///< Structure describing the material of the object
  /** A function computing an intersection, which returns the structure Hit */
  virtual Hit intersect(Ray ray) = 0;

  /** Function that returns the material struct of the object*/
  Material getMaterial() { return material; }
  /** Function that set the material
   @param material A structure describing the material of the object
  */
  void setMaterial(Material material) { this->material = material; }

  void setTransformation(glm::mat4 matrix)
  {
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
  Ray toLocal(const Ray &ray)
  {
    return Ray(
        vec3(inverseTransformationMatrix * glm::vec4(ray.origin, 1)),
        normalize(vec3(inverseTransformationMatrix * glm::vec4(ray.direction, 0))));
  }

  /**
   * Transform hit into global coords:
   * M * intersection and M^(-t) * normal to give it back in global coords
   * @param hit hit in local coords
   * @return hit in global coords
   */
  Hit toGlobal(const Hit &h, const vec3 &ray_origin)
  {
    return {
        .hit = h.hit,
        .normal = normalize(vec3(normalMatrix * glm::vec4(h.normal, 0))),
        .intersection = vec3(transformationMatrix * glm::vec4(h.intersection, 1)),
        .distance = glm::distance(h.intersection, ray_origin),
        .object = h.object,
        .uv = h.uv,
    };
  }
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object
{
private:
  float radius = 1; ///< Radius of the sphere
  glm::vec3 center; ///< Center of the sphere

public:
  /**
   The constructor of the sphere
   @param radius Radius of the sphere
   @param center Center of the sphere
   @param color Color of the sphere
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
  Hit intersect(Ray ray)
  {
    Ray local_ray = toLocal(ray);
    vec3 c = center - local_ray.origin;

    float cdotc = glm::dot(c, c);
    float cdotd = glm::dot(c, local_ray.direction);
    float D = 0;
    if (cdotc > cdotd * cdotd)
    {
      D = sqrt(cdotc - cdotd * cdotd);
    }

    if (D > radius)
    {
      return Hit();
    }

    float t1 = cdotd - sqrt(radius * radius - D * D);
    float t2 = cdotd + sqrt(radius * radius - D * D);
    float t = t1;
    if (t < 0)
      t = t2;

    if (t < 0)
    {
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
class Plane : public Object
{

private:
  glm::vec3 normal;
  glm::vec3 point;

public:
  Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal) {}
  Plane(glm::vec3 point, glm::vec3 normal, Material material)
      : point(point), normal(normal)
  {
    this->material = material;
  }
  Hit intersect(Ray ray)
  {
    Hit hit;
    hit.hit = false;
    float DdotN = glm::dot(ray.direction, normal);

    // always from behind, oriented from camera position to plane, not working
    // for oriented planes downwards
    if (DdotN < 0)
    {
      float PdotN = glm::dot(point - ray.origin, normal);
      float t = PdotN / DdotN;

      if (t > 0)
      {
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
  Plane *plane;

public:
  Cone(Material material)
  {
    this->material = material;
    plane = new Plane(glm::vec3(0, 1, 0), glm::vec3(0.0, 1, 0));
  }
  Hit intersect(Ray ray)
  {

    Hit hit;
    hit.hit = false;

    glm::vec3 d = inverseTransformationMatrix *
                  glm::vec4(ray.direction, 0.0); // implicit cast to vec3
    glm::vec3 o = inverseTransformationMatrix *
                  glm::vec4(ray.origin, 1.0); // implicit cast to vec3
    d = glm::normalize(d);

    // https://lousodrome.net/blog/light/2017/01/03/intersection-of-a-ray-and-a-cone/
    float a = d.x * d.x + d.z * d.z - d.y * d.y;
    float b = 2 * (d.x * o.x + d.z * o.z - d.y * o.y);
    float c = o.x * o.x + o.z * o.z - o.y * o.y;

    float delta = b * b - 4 * a * c;

    if (delta < 0)
    {
      return hit;
    }

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    float t = t1;
    hit.intersection = o + t * d;
    if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0)
    {
      t = t2;
      hit.intersection = o + t * d;
      if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0)
      {
        return hit;
      }
    };

    hit.normal =
        glm::vec3(hit.intersection.x, -hit.intersection.y, hit.intersection.z);
    hit.normal = glm::normalize(hit.normal);

    Ray new_ray(o, d);
    Hit hit_plane = plane->intersect(new_ray);

    // 3 cond: if point is inside circle base plane circumference (<= 1.0).
    // inside length(...) is the distance from intersection point to circle base
    // center
    if (hit_plane.hit && hit_plane.distance < t &&
        length(hit_plane.intersection - glm::vec3(0, 1, 0)) <= 1.0)
    {
      hit.intersection = hit_plane.intersection;
      hit.normal = hit_plane.normal;
    }

    hit.hit = true;
    hit.object = this;
    hit.intersection =
        transformationMatrix *
        glm::vec4(hit.intersection, 1.0); // implicit cast to vec3
    hit.normal =
        (normalMatrix * glm::vec4(hit.normal, 0.0)); // implicit cast to vec3
    hit.normal = glm::normalize(hit.normal);
    hit.distance = glm::length(hit.intersection - ray.origin);

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
  Light(glm::vec3 position) : position(position) { color = glm::vec3(1.0); }
  Light(glm::vec3 position, glm::vec3 color)
      : position(position), color(color) {}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001, 0.001, 0.001);
vector<Object *> objects; ///< A list of all objects in the scene

Hit get_closest_hit(Ray ray)
{
  Hit closest_hit;
  closest_hit.hit = false;
  closest_hit.distance = INFINITY;

  for (int k = 0; k < objects.size(); k++)
  {
    Hit h = objects[k]->intersect(ray);
    if (h.hit && h.distance < closest_hit.distance)
      closest_hit = h;
  }
  return closest_hit;
}

bool is_light_covered(vec3 intersection_point, vec3 intersection_normal,
                      vec3 light_direction, float light_distance)
{
  if (glm::dot(intersection_normal, light_direction) < 0)
  {
    return true;
  }

  Ray shadow =
      Ray(intersection_point + 0.001f * light_direction, light_direction);
  for (Object *object : objects)
  {
    Hit hit = object->intersect(shadow);
    if (hit.hit && hit.distance <= light_distance)
    {
      return true;
    }
  }
  return false;
}

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the
 viewer/camera
 @param material A material structure representing the material of the object
 @param recursion_depth A recursion depth of the raytracing
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv,
                     glm::vec3 view_direction, Material material,
                     int recursion_depth)
{
  glm::vec3 color(0.0);

  for (int light_num = 0; light_num < lights.size(); light_num++)
  {
    glm::vec3 light_direction =
        glm::normalize(lights[light_num]->position - point);

    float light_distance = glm::distance(lights[light_num]->position, point);
    if (is_light_covered(point, normal, light_direction, light_distance))
    {
      continue;
    }

    glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

    float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
    float VdotR =
        glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

    glm::vec3 diffuse_color = material.diffuse;
    if (material.texture)
    {
      diffuse_color = material.texture(uv);
    }

    glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
    glm::vec3 specular =
        material.specular * glm::vec3(pow(VdotR, material.shininess));

    // distance to the light
    float r = glm::distance(point, lights[light_num]->position);
    r = max(r, 0.1f);

    float shadow = 1.0f;
    Ray light_r = Ray(point + light_direction * 0.0001f, light_direction);
    for (int i = 0; i < objects.size(); i++)
    {
      Hit h = objects[i]->intersect(light_r);
      if (h.hit && h.distance < r)
      {
        shadow = 0.0;
        break;
      }
    }

    color += lights[light_num]->color * shadow * (diffuse + specular) / r / r;
  }

  // recursive raytracing
  if (recursion_depth > 0)
  {

    // Reflection
    glm::vec3 reflection(0.0f);
    if (material.reflection > 0)
    {
      color *= 1 - material.reflection;
      glm::vec3 reflection_direction = glm::reflect(-view_direction, normal);
      Ray reflection_ray =
          Ray(point + 0.001f * reflection_direction, reflection_direction);
      Hit closest_hit = get_closest_hit(reflection_ray);
      if (closest_hit.hit)
      {
        reflection = material.reflection *
                     PhongModel(
                         closest_hit.intersection, closest_hit.normal,
                         closest_hit.uv, glm::normalize(-reflection_direction),
                         closest_hit.object->getMaterial(), recursion_depth - 1);
      }
    }

    // Refraction
    glm::vec3 refraction(0.0f);
    if (material.refraction > 0)
    {
      color *= (1 - material.refraction);
      const bool is_entering = glm::dot(normal, -view_direction) < 0.0f;
      const float n1 =
          is_entering ? 1.0f : material.refraction_index; // assuming air
      const float n2 = is_entering ? material.refraction_index : 1.0f;
      const float beta = n1 / n2;
      glm::vec3 refraction_direction =
          glm::refract(-view_direction, is_entering ? normal : -normal, beta);
      Ray refraction_ray =
          Ray(point + 0.001f * refraction_direction, refraction_direction);

      Hit closest_hit = get_closest_hit(refraction_ray);
      if (closest_hit.hit)
      {
        refraction = material.refraction *
                     PhongModel(
                         closest_hit.intersection, closest_hit.normal,
                         closest_hit.uv, glm::normalize(-refraction_direction),
                         closest_hit.object->getMaterial(), recursion_depth - 1);
        // Fresnel effect
        float O1 = cos(glm::angle(normal, view_direction));
        float O2 = cos(glm::angle(-normal, refraction_direction));

        float R = 0.5 * (pow((n1 * O1 - n2 * O2) / (n1 * O1 + n2 * O2), 2) +
                         pow((n1 * O2 - n2 * O1) / (n1 * O2 + n2 * O1), 2));
        float T = 1 - R;

        reflection *= R;
        refraction *= T;
      }
    }
    color += reflection + refraction;
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
  Hit closest_hit = get_closest_hit(ray);
  glm::vec3 color(0.0);

  if (closest_hit.hit)
  {
    color = PhongModel(closest_hit.intersection, closest_hit.normal,
                       closest_hit.uv, glm::normalize(-ray.direction),
                       closest_hit.object->getMaterial(), 3);
  }
  // clamp the final color to [0,1]
  return glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
}
/**
 Function defining the scene
 */
void sceneDefinition(float ya, float yb)
{

  Material green_diffuse;
  green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
  green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

  Material red_specular;
  red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
  red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
  red_specular.specular = glm::vec3(0.5);
  red_specular.shininess = 10.0;

  Material blue_specular;
  blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
  blue_specular.diffuse = glm::vec3(0.0f, 0.0f, 0.0f);
  blue_specular.specular = glm::vec3(0.6);
  blue_specular.shininess = 100.0;
  blue_specular.reflection = 1.0f;

  Material refractive;
  //  refractive.ambient = glm::vec3(0.5f);
  //  refractive.diffuse = glm::vec3(0.5f);
  refractive.specular = glm::vec3(1.0f);
  refractive.shininess = 10.0f;
  refractive.reflection = 1.0f;
  refractive.refraction = 1.0f;
  refractive.refraction_index = 2.0f;

  Sphere *s0 = new Sphere(refractive);

  s0->setTransformation(
      translate(vec3(-3.0f, (yb + 1), 8)) *
      glm::scale(vec3(2.0f)));
  objects.push_back(s0);

  Sphere *s1 = new Sphere(blue_specular);
  if (ya > -2 && ya < -1)
  {
    s1->setTransformation(glm::translate(glm::vec3(1, ya, 8)) *
                          glm::scale(glm::vec3(1.0f, 1 / (ya * (-1)), 1.0f)));
  }
  else
  {
    s1->setTransformation(glm::translate(glm::vec3(1, ya, 8)) *
                          glm::scale(glm::vec3(1.0f, 1.0f, 1.0f)));
  }

  objects.push_back(s1);

  Sphere *s2 = new Sphere(red_specular);
  s2->setTransformation(
      translate(vec3(-0.7, -1, 3.5)) *
      glm::scale(glm::mat4(1), vec3(0.5)));
  // objects.push_back(new Sphere(1.0, glm::vec3(3,-2,6), green_diffuse));

  // Textured sphere

  Material textured;
  textured.texture = &rainbowTexture;

  Sphere *s3 = new Sphere(textured);
  s3->setTransformation(
      translate(glm::vec3(-6, 4, 23)) *
      scale(glm::vec3(7)));

  objects.push_back(s3);

  // Planes
  Material red_diffuse;
  red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
  red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);

  Material blue_diffuse;
  blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
  blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
  objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0.0, 1, 0)));
  objects.push_back(
      new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), green_diffuse));
  objects.push_back(
      new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));
  objects.push_back(
      new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));
  objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));
  objects.push_back(new Plane(glm::vec3(0, 1, -0.01), glm::vec3(0.0, 0.0, 1.0),
                              green_diffuse));

  // Cone

  Material yellow_specular;
  yellow_specular.ambient = glm::vec3(0.1f, 0.10f, 0.0f);
  yellow_specular.diffuse = glm::vec3(0.4f, 0.4f, 0.0f);
  yellow_specular.specular = glm::vec3(1.0);
  yellow_specular.shininess = 100.0;

  Cone *cone = new Cone(yellow_specular);
  glm::mat4 translationMatrix = glm::translate(glm::vec3(5, 9, 14));
  glm::mat4 scalingMatrix = glm::scale(glm::vec3(3.0f, 12.0f, 3.0f));
  glm::mat4 rotationMatrix =
      glm::rotate(glm::radians(180.0f), glm::vec3(1, 0, 0));
  cone->setTransformation(translationMatrix * scalingMatrix * rotationMatrix);
  objects.push_back(cone);

  Cone *cone2 = new Cone(green_diffuse);
  translationMatrix = glm::translate(glm::vec3(6, -3, 7));
  scalingMatrix = glm::scale(glm::vec3(1.0f, 3.0f, 1.0f));
  rotationMatrix = glm::rotate(glm::atan(3.0f), glm::vec3(0, 0, 1));
  cone2->setTransformation(translationMatrix * rotationMatrix * scalingMatrix);
  objects.push_back(cone2);

  lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
  lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
  lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
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
  return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)),
                    glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char *argv[])
{

  clock_t t = clock(); // variable for keeping the time of the rendering

  int width = 1024; // width of the image
  int height = 768; // height of the image
  float fov = 90;   // field of view

  if (argc == 4)
  {
    float ya = stof(argv[2]);
    float yb = stof(argv[3]);

    sceneDefinition(ya, yb); // Let's define a scene
  }
  else
  {
    sceneDefinition(-2, 0); // Let's define a scene
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
  cout << "It took " << ((float)t) / CLOCKS_PER_SEC
       << " seconds to render the image." << endl;
  cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t)
       << " frames per second." << endl;

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
