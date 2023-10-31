/**
 * @file main.cpp
 */

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "glm/glm.hpp"

#include "Image.h"
#include "Material.h"

// #define PI atan(1) * 4

using glm::distance;
using glm::normalize;
using glm::reflect;
using glm::vec3;

// #region RAY
/**
 * Class representing a single ray.
 */
class Ray {
public:
  vec3 origin;
  vec3 direction;

  /**
   * Constructor of the ray
   * @param origin Origin of the ray
   * @param direction Direction of the ray
   */
  Ray(vec3 origin, vec3 direction) : origin(origin), direction(direction) {}
};
// #endregion

// #region HIT
class Object;

/**
 * Structure representing the event of hitting an object
 */
struct Hit {
  bool hit = false; // state if intersection with an object
  vec3 normal; // Normal vector of the intersected object at the intersection
               // point
  vec3 intersection;         // Point of Intersection
  float distance = INFINITY; // Distance from the origin of the ray to the
                             // intersection point
  Object *object = 0;        // pointer to the intersected object
  vec2 uv; // Coordinates for computing the texture (texture coordinates)
};
// #endregion

// #region OBJECT
/**
 * General class for the object
 */
class Object {
public:
  vec3 color = vec3(0.0f);
  Material material; // Structure describing the material of the object

  /** A function computing an intersection, which returns the structure Hit */
  virtual Hit intersect(Ray ray) = 0;

  /** Function that returns the material struct of the object*/
  Material getMaterial() const { return material; }

  /**
   * Function that set the material
   * @param material A structure describing the material of the object
   */
  void setMaterial(Material m) { this->material = m; }
};
// #endregion

// #region SPHERE
/**
 * Implementation of the class Object for sphere shape.
 */
class Sphere : public Object {
private:
  float radius;
  vec3 center;

public:
  /**
   * Constructor of the sphere
   * @param radius Radius of the sphere
   * @param center Center of the sphere
   * @param color Color of the sphere
   */
  Sphere(float radius, vec3 center, vec3 color)
      : radius(radius), center(center) {
    this->color = color;
  }

  Sphere(float radius, vec3 center, Material material)
      : radius(radius), center(center) {
    this->material = material;
  }

  /**
   * Implementation of the intersection function
   */
  Hit intersect(Ray ray) {
    vec3 c = center - ray.origin;
    float cdotc = glm::dot(c, c);
    float cdotd = glm::dot(c, ray.direction);
    Hit hit;
    float D = 0;
    if (cdotc > cdotd * cdotd) {
      D = sqrt(cdotc - cdotd * cdotd);
    }

    if (D > radius) {
      hit.hit = false;
      return hit;
    }

    hit.hit = true;
    float t1 = cdotd - sqrt(radius * radius - D * D);
    float t2 = cdotd + sqrt(radius * radius - D * D);

    float t = t1;
    if (t < 0)
      t = t2;

    if (t < 0) {
      hit.hit = false;
      return hit;
    }

    hit.intersection = ray.origin + t * ray.direction;
    hit.normal = normalize(hit.intersection - center);
    hit.distance = distance(ray.origin, hit.intersection);
    hit.object = this;

    // Ex2: computing texture coordinates for the sphere.
    hit.uv.s = (asin(hit.normal.y) + M_PI / 2) / M_PI;
    hit.uv.t = (atan2(hit.normal.z, hit.normal.x) + M_PI) / (2 * M_PI);
    return hit;
  }
};
// #endregion

// #region PLANE & TRIANGLE
class Plane : public Object {
private:
  vec3 point = vec3(0);
  vec3 normal = vec3(0, 0, 1);

public:
  Plane(vec3 point, vec3 normal) : point(point), normal(normal) {}
  Plane(vec3 point, vec3 normal, Material material)
      : point(point), normal(normal) {
    this->material = material;
  }

  /**
   * Ex1: Plane-ray intersection
   */
  Hit intersect(Ray ray) {
    float ray_dot_n = dot(point - ray.origin, normal);
    float d_dot_N = dot(ray.direction, normal);

    // cos b of angle parallel to plane
    if (d_dot_N == 0)
      return Hit();

    float t = ray_dot_n / d_dot_N;

    // intersection behind ray origin
    if (t < 0)
      return Hit();

    Hit hit{
        .hit = true,
        .normal = normalize(normal),
        .intersection = ray.origin + t * ray.direction,
        .distance = distance(ray.origin, hit.intersection),
        .object = this,
        .uv = {0, 0},
    };
    return hit;
  }
};

class Triangle : public Object {
private:
  vec3 point_a, point_b, point_c;
  vec3 normal_a = vec3(0), normal_b = vec3(0), normal_c = vec3(0);
  bool with_normal = false;

public:
  Triangle(vec3 pa, vec3 pb, vec3 pc, Material mat)
      : point_a(pa), point_b(pb), point_c(pc) {
    this->setMaterial(mat);
  }

  Triangle(vec3 pa, vec3 pb, vec3 pc, vec3 na, vec3 nb, vec3 nc, Material mat)
      : point_a(pa), point_b(pb), point_c(pc), normal_a(na), normal_b(nb),
        normal_c(nc) {
    this->setMaterial(mat);
    with_normal = true;
  }

  Hit intersect(Ray ray) override {
    vec3 p1 = point_a, p2 = point_b, p3 = point_c;
    vec3 edge1 = p2 - p1;
    vec3 edge2 = p3 - p1;
    vec3 perpendicular = cross(edge1, edge2);
    Hit plane_hit = Plane(p1, normalize(perpendicular)).intersect(ray);

    // t check done inside the plane
    if (!plane_hit.hit)
      return {};

    // check if p inside triangle
    vec3 p = plane_hit.intersection;
    float W = 0.5f * length(perpendicular);
    auto sign = [](float x) { return 0 <= x ? 1.0f : -1.0f; };

    vec3 perpendicular1 = cross(p2 - p, p3 - p);
    float sign1 = sign(dot(perpendicular1, perpendicular));
    float w1 = 0.5f * length(perpendicular1) * sign1;

    vec3 perpendicular2 = cross(p3 - p, p1 - p);
    float sign2 = sign(dot(perpendicular2, perpendicular));
    float w2 = 0.5f * length(perpendicular2) * sign2;

    vec3 perpendicular3 = cross(p1 - p, p2 - p);
    float sign3 = sign(dot(perpendicular3, perpendicular));
    float w3 = 0.5f * length(perpendicular3) * sign3;

    if (1.0f != sign1 || sign1 != sign2 || sign2 != sign3)
      return {};

    vec3 normal;
    if (with_normal)
      normal = normalize(w1 * normal_a + w2 * normal_b + w3 * normal_c);
    else
      normal = normalize(perpendicular);

    Hit hit = {
        .hit = true,
        .normal = normal,
        .intersection = p,
        .distance = plane_hit.distance,
        .object = this,
    };
    return hit;
  }
};

class Mesh : public Object {
public:
  vector<Triangle> triangles;
  string fname;
  string object_name;
  // no smoothing

public:
  Mesh(string fname_, vec3 translation, Material material) : fname(fname_) {
    setMaterial(material);
    load_mesh(translation);
  }

  void load_mesh(vec3 translation) {
    ifstream objFile(this->fname);
    string line;
    vector<vec3> v_positions;
    vector<vec3> vn_positions;

    while (getline(objFile, line)) {
      istringstream stream(line);
      string token;
      stream >> token;

      if (token == "v") {
        float x, y, z;
        stream >> x >> y >> z;
        v_positions.push_back(vec3(x, y, z) + translation);
      } else if (token == "vn") {
        // normals might not be unit vectors
        float x, y, z;
        stream >> x >> y >> z;
        vn_positions.push_back(normalize(vec3(x, y, z)));
      } else if (token == "f") {
        string i, j, k;
        stream >> i >> j >> k;

        bool has_normals = i.find("//") != string::npos;
        if (has_normals) {
          int i_v = stoi(i.substr(0, i.find("//")));
          int i_vn = stoi(i.substr(i.find("//") + 2));

          int j_v = stoi(j.substr(0, j.find("//")));
          int j_vn = stoi(j.substr(j.find("//") + 2));

          int k_v = stoi(k.substr(0, k.find("//")));
          int k_vn = stoi(k.substr(k.find("//") + 2));

          this->triangles.push_back(Triangle(
              v_positions[i_v - 1], v_positions[j_v - 1], v_positions[k_v - 1],
              vn_positions[i_vn - 1], vn_positions[j_vn - 1],
              vn_positions[k_vn - 1], this->material));
        } else {
          this->triangles.push_back(
              Triangle(v_positions[stoi(i) - 1], v_positions[stoi(j) - 1],
                       v_positions[stoi(k) - 1], this->material));
        }
      } else if (token == "o") {
        stream >> this->object_name;
      }
      // add smoothing checks
    }
    objFile.close();
  }

  Hit intersect(Ray ray) override {
    Hit closest_hit;
    for (Triangle t : this->triangles) {
      Hit h = t.intersect(ray);
      if (h.hit && h.distance < closest_hit.distance)
        closest_hit = h;
    }
    closest_hit.object = this;
    return closest_hit;
  }
};

// #region LIGHT
/**
 * Light class
 */
class Light {
public:
  vec3 position;
  vec3 color;

  Light(vec3 position, vec3 color) : position(position), color(color) {}

  explicit Light(vec3 position) { Light(position, vec3(1.0f)); }
};
// #endregion

// GLOBAL VARIABLES
vector<Light *> lights; // list of lights in the scene
vec3 ambient_light(0.5f);
vector<Object *> objects; // list of all objects in the scene

// #region PHONG-MODEL
/**
 * Function for computing color of an object according to the Phong Model
 * @param point A point belonging to the object for which the color is computed
 * @param normal A normal vector the the point
 * @param uv Texture coordinates
 * @param view_direction A normalized direction from the point to the
 * viewer/camera
 * @param material A material structure representing the material of the object
 */
vec3 PhongModel(vec3 point, vec3 normal, vec2 uv, vec3 view_direction,
                Material material) {
  vec3 color(0.0f);

  vec3 I_diffuse = vec3(0), I_specular = vec3(0),
       // ambient illumination
      I_ambient = material.ambient * ambient_light;

  for (Light *&l : lights) {
    /* Ex3: Modify the code by adding attenuation of the light due to distance
     * from the intersection point to the light source
     */
    float r = distance(point, l->position);
    float attenuation = 1 / pow(max(r, 0.5f), 2);

    // DIFFUSE REFLECTION:
    // from diffuse reflection point, cos of angle of surface normal n and
    // direction from point to light source
    vec3 l_direction = normalize(l->position - point);
    float cos_phi =
        glm::clamp(dot(normal, l_direction), 0.0f, 1.0f); // already normalized

    /* Ex2: Modify the code by adding texturing,
     * i.e. diffuse color should be computed using one of the texture functions
     * according to the texture coordinates stored in the uv variable.
     * Make sure that the code works also for objects that should not have
     * texture.
     */

    I_diffuse += material.diffuse * cos_phi * l->color * attenuation *
                 (material.texture ? material.texture(uv) : vec3(1.0f));

    // SPECULAR HIGHLIGHT:
    vec3 r_direction = reflect(-l_direction, normal);
    float k = max(material.shininess, 1.0f);
    float cos_alpha = glm::clamp(dot(view_direction, r_direction), 0.0f, 1.0f);
    I_specular +=
        material.specular * pow(cos_alpha, k) * l->color * attenuation;
  }

  color = I_ambient + I_diffuse + I_specular;
  return clamp(color, vec3(0), vec3(1));
}
// #endregion

// #region TRACE-RAY
/**
 * Function computing a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 trace_ray(Ray ray) {
  Hit closest_hit;

  for (Object *&o : objects) {
    Hit hit = o->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance)
      closest_hit = hit;
  }

  if (!closest_hit.hit) {
    return vec3(0);
  }
  Object *obj = closest_hit.object;
//  cout << "obj " << obj << endl;
  Mesh *mesh = dynamic_cast<Mesh *>(obj);
  if (mesh) {
//  cout << mesh->object_name << " " << mesh << endl;
//  cout << "Material: " << obj->getMaterial().diffuse.x << " "
//       << obj->getMaterial().diffuse.y << " " << obj->getMaterial().diffuse.z
//       << endl;

  }
  return PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv,
               normalize(-ray.direction), closest_hit.object->getMaterial());
}
// #endregion

// #region SCENE
/**
 * Function defining the scene
 */
void sceneDefinition() {
  vec3 color_red{1.5f, 0.3f, 0.3f};
  vec3 color_blue{0.4f, 0.4f, 1.5f};
  vec3 color_green{0.5f, 1.5f, 0.5f};
  vec3 color_white{1.5f, 1.5f, 1.5f};

  Material red_specular{
      .ambient = {0.01f, 0.03f, 0.03f},
      .diffuse = color_red,
      .specular = vec3(0.5f),
      .shininess = 10.0f,
  };

  Material blue_shiny{
      .ambient = {0.07f, 0.07f, 0.1f},
      .diffuse = color_blue,
      .specular = vec3(0.6f),
      .shininess = 100.0f,
  };

  Material green_diffuse{
      .ambient = {0.07f, 0.09f, 0.07f},
      .diffuse = color_green,
  };

  Material white_plain{
      .ambient = {1.0f, 1.0f, 1.0f},
      .diffuse = color_white,
      .specular = vec3(1.6f),
      .shininess = 0.0f,
  };

  Material material_rainbow{
      .texture = &rainbowTexture,
  };

  // TEST
  objects.push_back(new Mesh("./meshes/armadillo_with_normals.obj",
                             {-4, -3, 10},
                             {
                                 .ambient = {0.7f, 0.7f, 0.1f},
                                 .diffuse = {0.3, 0.3, 1.0f},
                                 .specular = {0.2f, 0.2f, 0.5f},
                                 .shininess = 40.0f,
                             }));
  objects.push_back(new Mesh("./meshes/lucy_with_normals.obj", {4, -3, 10},
                             {
                                 .ambient = {0.2f, 0.5f, 0.2f},
                                 .diffuse = {0.4f, 1.0f, 0.4f},
                                 .specular = {0.2f, 0.4f, 0.1f},
                                 .shininess = 30.0f,
                             }));
  objects.push_back(new Mesh("./meshes/bunny_with_normals.obj", {0, -3, 8},
                             {
                                 .ambient = {0.7f, 0.1f, 0.1f},
                                 .diffuse = {1.0f, 0.3f, 0.3f},
                                 .specular = {0.4f, 0.1f, 0.1f},
                                 .shininess = 50.0f,
                             }));
  // objects.push_back(new Mesh("./meshes/bunny.obj",
  //                             {0, -2, 8},
  //                             {.diffuse = {0.25f, 0.25f, 0.5f},}));
  // objects.push_back(new Mesh("./meshes/"))

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
    objects.push_back(new Plane({0.0f, 0.0f, -0.01f}, {0.0f, 0.0f, 1.0f},
                                Material())); // back
    objects.push_back(new Plane({15.0f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f},
                                {
                                    .diffuse = {0.6f, 0.6f, 1.0f},
                                })); // right
    objects.push_back(new Plane({0.0f, 0.0f, 30.0f}, {0.0f, 0.0f, -1.0f},
                                {
                                    .diffuse = {0.5f, 1.0f, 0.5f},
                                })); // front
    objects.push_back(new Plane({-15.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f},
                                {
                                    .diffuse = {0.6f, 0.4f, 0.4f},
                                })); // left
    objects.push_back(new Plane({0.0f, 27.0f, 0.0f}, {0.0f, -1.0f, 0.0f},
                                Material())); // top
    objects.push_back(new Plane({0.0f, -3.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                Material())); // bottom

  // lights
  lights.push_back(new Light({0.0f, 26.0f, 2.0f}, vec3(100.0f)));
  lights.push_back(new Light({0.0f, 5.0f, 1.0f}, vec3(40.0f)));
  lights.push_back(new Light({0.0f, 1.0f, 3.0f}, vec3(25.0f)));
}
// #endregion

// #region TONE-MAPPING
/**
 * Performing tone mapping and gamma correction of intensity computed using the
 * raytracer EX3 value of gamma from:
 * https://www.rtings.com/laptop/reviews/dell/precision-5560-2021#test_5194
 * @param intensity Input intensity
 * @return Tone mapped intensity in range [0,1]
 */
vec3 toneMapping(vec3 intensity) {
  float alpha = 0.8f, beta = 1.2f, gamma = 2.13f;
  vec3 I_tone_mapped = alpha * pow(intensity, vec3(beta));
  vec3 I_gamma_corrected = min(pow(I_tone_mapped, vec3(1 / gamma)), 1.0f);
  return clamp(I_gamma_corrected, vec3(0.0f), vec3(1.0f));
}
// #endregion

int main(int argc, char const *argv[]) {
  clock_t t = clock(); // keeping the time of the rendering
                       // Default
                         int width = 1024;
                         int height = 768;
                       // Debug
//  int width = 320;
//  int height = 210;
  // Final
//   int width = 2048;
//   int height = 1536;
  float fov = 90; // field of view
  sceneDefinition();
  Image image(width, height); // Create an image where we will store the result
  float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
  float X = -s * width / 2;
  float Y = s * height / 2;

  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++) {
      float dx = X + i * s + s / 2;
      float dy = Y - j * s - s / 2;
      vec3 direction(dx, dy, 1);
      direction = normalize(direction);
      Ray ray(vec3(0), direction);
      image.setPixel(i, j, toneMapping(trace_ray(ray)));
    }

  t = clock() - t;
  cout << "It took " << (float)t / CLOCKS_PER_SEC
       << " seconds to render the image." << endl;
  cout << "I could render at " << (float)CLOCKS_PER_SEC / (float)t
       << " frames per second." << endl;

  // Writing the final results of the rendering
  if (argc == 2)
    image.writeImage(argv[1]);
  else
    image.writeImage("./result.ppm");

  return 0;
}
