/**
 * @file main.cpp
 */
#define DEBUG 0

#include "glm/glm.hpp"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <vector>

#include "Image.h"
#include "Material.h"

#include <climits>

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
  virtual Hit intersect(Ray &ray) = 0;

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
  Hit intersect(Ray &ray) override {
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
  Hit intersect(Ray &ray) override {
    float ray_dot_n = dot(point - ray.origin, normal);
    float d_dot_N = dot(ray.direction, normal);

    // cos b of angle parallel to plane
    if (d_dot_N == 0)
      return {};

    float t = ray_dot_n / d_dot_N;

    // intersection behind ray origin
    if (t < 0)
      return {};

    Hit hit;
    hit.hit = true;
    hit.normal = normalize(normal);
    hit.intersection = ray.origin + t * ray.direction;
    hit.distance = distance(ray.origin, hit.intersection);
    hit.object = this;
    hit.uv = {0, 0};

    return hit;
  }
};

class Triangle : public Object {

public:
  vec3 point_a, point_b, point_c;
  vec3 normal_a = vec3(0), normal_b = vec3(0), normal_c = vec3(0);
  bool with_normal = false;
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

  Hit intersect(Ray &ray) override {
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

    Hit hit;
    hit.hit = true;
    hit.normal = normal;
    hit.intersection = p;
    hit.distance = plane_hit.distance;
    hit.object = this;

    return hit;
  }
};

bool bboxIntersect(Ray &ray, float mi_x, float ma_x, float mi_y, float ma_y,
                   float mi_z, float ma_z) {
  float tmin = (mi_x - ray.origin.x) / ray.direction.x;
  float tmax = (ma_x - ray.origin.x) / ray.direction.x;

  if (tmin > tmax)
    std::swap(tmin, tmax);

  float tymin = (mi_y - ray.origin.y) / ray.direction.y;
  float tymax = (ma_y - ray.origin.y) / ray.direction.y;

  if (tymin > tymax)
    std::swap(tymin, tymax);

  if ((tmin > tymax) || (tymin > tmax))
    return false;

  if (tymin > tmin)
    tmin = tymin;

  if (tymax < tmax)
    tmax = tymax;

  float tzmin = (mi_z - ray.origin.z) / ray.direction.z;
  float tzmax = (ma_z - ray.origin.z) / ray.direction.z;

  if (tzmin > tzmax)
    std::swap(tzmin, tzmax);

  if ((tmin > tzmax) || (tzmin > tmax))
    return false;

  return true;
}

struct AABB {
  glm::vec3 min;
  glm::vec3 max;
};

struct BVHNode {
  AABB bounds;
  vector<Triangle> triangles;
  BVHNode *left;
  BVHNode *right;
};

AABB calculateBoundingBox(const std::vector<Triangle> &triangles) {
  AABB box;
  box.min = vec3(INT_MAX);
  box.max = vec3(INT_MIN);

  for (const Triangle &triangle : triangles) {
    box.min.x = std::min(
        box.min.x, std::min(triangle.point_a.x,
                            std::min(triangle.point_b.x, triangle.point_c.x)));
    box.min.y = std::min(
        box.min.y, std::min(triangle.point_a.y,
                            std::min(triangle.point_b.y, triangle.point_c.y)));
    box.min.z = std::min(
        box.min.z, std::min(triangle.point_a.z,
                            std::min(triangle.point_b.z, triangle.point_c.z)));

    box.max.x = std::max(
        box.max.x, std::max(triangle.point_a.x,
                            std::max(triangle.point_b.x, triangle.point_c.x)));
    box.max.y = std::max(
        box.max.y, std::max(triangle.point_a.y,
                            std::max(triangle.point_b.y, triangle.point_c.y)));
    box.max.z = std::max(
        box.max.z, std::max(triangle.point_a.z,
                            std::max(triangle.point_b.z, triangle.point_c.z)));
  }

  return box;
}

std::pair<std::vector<Triangle>, std::vector<Triangle>>
splitTrianglesSpace(std::vector<Triangle> &triangles, int axis) {
  std::vector<Triangle> left;
  std::vector<Triangle> right;

  float mid = 0;
  for (const Triangle &triangle : triangles) {
    mid += triangle.point_a[axis] + triangle.point_b[axis] +
           triangle.point_c[axis];
  }
  mid /= triangles.size() * 3;

  for (const Triangle &triangle : triangles) {
    if (triangle.point_a[axis] < mid || triangle.point_b[axis] < mid ||
        triangle.point_c[axis] < mid) {
      left.emplace_back(triangle);
    } else {
      right.emplace_back(triangle);
    }
  }

  return {left, right};
}

// Function to build the BVH recursively
BVHNode *buildBVH(std::vector<Triangle> &triangles, int axis = 0) {
  auto *node = new BVHNode;
  node->bounds = calculateBoundingBox(triangles);

  if (triangles.size() <= 20) {
    node->left = nullptr;
    node->right = nullptr;
    node->triangles = triangles;
  } else {
    std::pair<std::vector<Triangle>, std::vector<Triangle>> split =
        splitTrianglesSpace(triangles, axis);
    std::vector<Triangle> &left = split.first;
    std::vector<Triangle> &right = split.second;

    node->left = buildBVH(left, (axis + 1) % 3);
    node->right = buildBVH(right, (axis + 1) % 3);
  }

  return node;
}

vector<Triangle> intersectBVH(BVHNode *node, Ray &ray) {

  if (node->left == nullptr && node->right == nullptr) {
    return node->triangles;
  }

  // If right and left, return both
  if (bboxIntersect(ray, node->left->bounds.min.x, node->left->bounds.max.x,
                    node->left->bounds.min.y, node->left->bounds.max.y,
                    node->left->bounds.min.z, node->left->bounds.max.z)

      &&
      bboxIntersect(ray, node->right->bounds.min.x, node->right->bounds.max.x,
                    node->right->bounds.min.y, node->right->bounds.max.y,
                    node->right->bounds.min.z, node->right->bounds.max.z)) {
    vector<Triangle> left = intersectBVH(node->left, ray);
    vector<Triangle> right = intersectBVH(node->right, ray);
    left.insert(left.end(), right.begin(), right.end());
    return left;
  }
  // if left, return left
  if (bboxIntersect(ray, node->left->bounds.min.x, node->left->bounds.max.x,
                    node->left->bounds.min.y, node->left->bounds.max.y,
                    node->left->bounds.min.z, node->left->bounds.max.z)) {
    return intersectBVH(node->left, ray);
  }
  // if right, return right
  if (bboxIntersect(ray, node->right->bounds.min.x, node->right->bounds.max.x,
                    node->right->bounds.min.y, node->right->bounds.max.y,
                    node->right->bounds.min.z, node->right->bounds.max.z)) {
    return intersectBVH(node->right, ray);
  }

#if DEBUG
  cout << "EMPTY, is this normal?" << endl;
  cout << node->triangles.size() << endl;
#endif
  return node->triangles;
}
class Mesh : public Object {
public:
  vector<Triangle> triangles;
  string fname;
  string object_name;

  float min_x = INT_MAX;
  float max_x = INT_MIN;

  float min_y = INT_MAX;
  float max_y = INT_MIN;

  float min_z = INT_MAX;
  float max_z = INT_MIN;

  bool is_bvh;
  BVHNode *bVHNode;

public:
  Mesh(string fname_, bool is_bvh_, vec3 translation, Material material)
      : fname(fname_), is_bvh(is_bvh_) {
    setMaterial(material);
    load_mesh(translation);
#if DEBUG
    cout << min_x << " " << min_y << " " << min_z << endl;
    cout << max_x << " " << max_y << " " << max_z << endl;
    cout << "is BVH = " << is_bvh << endl;
#endif

    if (this->is_bvh) {
#if DEBUG
      cout << "BUILDING BVH" << endl;
#endif
      bVHNode = buildBVH(this->triangles);
#if DEBUG
      cout << "BUILT BVH" << endl;
#endif
    }
  }

  void load_mesh(vec3 &translation) {
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
        vec3 pos = vec3(x, y, z) + translation;
        this->min_x = min(this->min_x, pos.x);
        this->max_x = max(this->max_x, pos.x);

        this->min_y = min(this->min_y, pos.y);
        this->max_y = max(this->max_y, pos.y);

        this->min_z = min(this->min_z, pos.z);
        this->max_z = max(this->max_z, pos.z);
        v_positions.push_back(pos);
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

          this->triangles.emplace_back(
              v_positions[i_v - 1], v_positions[j_v - 1], v_positions[k_v - 1],
              vn_positions[i_vn - 1], vn_positions[j_vn - 1],
              vn_positions[k_vn - 1], this->material);
        } else {
          this->triangles.emplace_back(
              v_positions[stoi(i) - 1], v_positions[stoi(j) - 1],
              v_positions[stoi(k) - 1], this->material);
        }
      }
    }
    objFile.close();

    cout << "Number of Triangles for " << this->fname << " = "
         << this->triangles.size() << endl;
  }

  Hit intersect(Ray &ray) override {
    Hit closest_hit;

    if (!bboxIntersect(ray, min_x, max_x, min_y, max_y, min_z, max_z))
      return closest_hit;

    if (this->is_bvh) {
      for (Triangle &t : intersectBVH(this->bVHNode, ray)) {
        Hit h = t.intersect(ray);
        if (h.hit && h.distance < closest_hit.distance)
          closest_hit = h;
      }
      closest_hit.object = this;
      return closest_hit;
    }

    for (Triangle &t : this->triangles) {
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
  explicit Light(vec3 position) : Light(position, vec3(1.0f)) {}
};
// #endregion

// GLOBAL VARIABLES
vector<Light *> lights; // list of lights in the scene
vec3 ambient_light(0.5f);
vector<Object *> objects; // list of all objects in the scene

// #region PHONG-MODEL
bool is_light_covered(vec3 intersection_point, vec3 intersection_normal,
                      vec3 light_direction, float light_distance) {
  if (glm::dot(intersection_normal, light_direction) < 0) {
    return true;
  }

  Ray shadow =
      Ray(intersection_point + 0.001f * light_direction, light_direction);
  for (Object *object : objects) {
    Hit hit = object->intersect(shadow);
    if (hit.hit && hit.distance <= light_distance) {
      return true;
    }
  }
  return false;
}

/**
 * Function for computing color of an object according to the Phong Model
 * @param point A point belonging to the object for which the color is computed
 * @param normal A normal vector the the point
 * @param uv Texture coordinates
 * @param view_direction A normalized direction from the point to the
 * viewer/camera
 * @param material A material structure representing the material of the object
 */
vec3 PhongModel(vec3 &point, vec3 &normal, vec2 &uv, vec3 &&view_direction,
                Material &&material) {
  glm::vec3 color(0.0);

  for (auto &light : lights) {
    glm::vec3 light_direction = glm::normalize(light->position - point);

    float light_distance = glm::distance(light->position, point);
    if (is_light_covered(point, normal, light_direction, light_distance)) {
      continue;
    }

    glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

    float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
    float VdotR =
        glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

    glm::vec3 diffuse_color = material.diffuse;
    if (material.texture) {
      diffuse_color = material.texture(uv);
    }

    glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
    //    glm::vec3 specular =
    //        material.specular * glm::vec3(pow(VdotR, material.shininess));

    // distance to the light
    float r = glm::distance(point, light->position);
    r = max(r, 0.1f);

    float shadow = 1.0f;
    Ray light_r = Ray(point + light_direction * 0.0001f, light_direction);
    for (auto &object : objects) {
      Hit h = object->intersect(light_r);
      if (h.hit && h.distance < r) {
        shadow = 0.0;
        break;
      }
    }

    color += light->color * shadow * (diffuse /*+ specular*/) / r / r;
  }
  color += ambient_light * material.ambient;
  return color;
}
// #endregion

// #region TRACE-RAY
/**
 * Function computing a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 trace_ray(Ray &ray) {
  Hit closest_hit;

  for (Object *&o : objects) {
    Hit hit = o->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance)
      closest_hit = hit;
  }

  if (!closest_hit.hit) {
    return vec3(0);
  }
  return PhongModel(closest_hit.intersection, closest_hit.normal,
                    closest_hit.uv, normalize(-ray.direction),
                    closest_hit.object->getMaterial());
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

  Material red_specular;
  red_specular.ambient = {0.01f, 0.03f, 0.03f};
  red_specular.diffuse = color_red;
  red_specular.specular = vec3(0.5f);
  red_specular.shininess = 10.0f;

  Material blue_shiny;
  blue_shiny.ambient = {0.07f, 0.07f, 0.1f};
  blue_shiny.diffuse = color_blue;
  blue_shiny.specular = vec3(0.6f);
  blue_shiny.shininess = 100.0f;

  Material green_diffuse;
  green_diffuse.ambient = {0.07f, 0.09f, 0.07f};
  green_diffuse.diffuse = color_green;

  Material white_plain;
  white_plain.ambient = {0.5f, 0.5f, 0.5f};
  white_plain.diffuse = color_white;
  white_plain.specular = vec3(1.6f);
  white_plain.shininess = 0.0f;

  // Objects
//   objects.push_back(new Mesh("./meshes/armadillo_small.obj", true, {-4, -3, 10}, white_plain));
//   objects.push_back(new Mesh("./meshes/lucy_small.obj", true, {4, -3, 10}, white_plain));
//   objects.push_back(new Mesh("./meshes/bunny_small.obj", true, {0, -3, 8}, white_plain));

  objects.push_back(
      new Mesh("./meshes/armadillo.obj", true, {-4, -3, 10}, white_plain));
  objects.push_back(
      new Mesh("./meshes/lucy.obj", true, {4, -3, 10}, white_plain));
  objects.push_back(
      new Mesh("./meshes/bunny.obj", true, {0, -3, 8}, white_plain));

  // planes
  objects.push_back(
      new Plane({0.0f, 0.0f, -0.01f}, {0.0f, 0.0f, 1.0f}, white_plain)); // back
  objects.push_back(new Plane({15.0f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f},
                              white_plain)); // right
  objects.push_back(new Plane({0.0f, 0.0f, 30.0f}, {0.0f, 0.0f, -1.0f},
                              white_plain)); // front
  objects.push_back(
      new Plane({-15.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, white_plain)); // left
  objects.push_back(
      new Plane({0.0f, 27.0f, 0.0f}, {0.0f, -1.0f, 0.0f}, white_plain)); // top
  objects.push_back(new Plane({0.0f, -3.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                              white_plain)); // bottom

  // lights
  lights.push_back(new Light({0.0f, 5.0f, 15.0f}, vec3(30.0f)));
  lights.push_back(new Light({0.0f, 5.0f, 1.0f}, vec3(15.0f)));
  lights.push_back(new Light({0.0f, 1.0f, 3.0f}, vec3(7.0f)));
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
vec3 toneMapping(vec3 &&intensity) {
  float alpha = 0.8f, beta = 1.2f, gamma = 2.13f;
  vec3 I_tone_mapped = alpha * pow(intensity, vec3(beta));
  vec3 I_gamma_corrected = min(pow(I_tone_mapped, vec3(1 / gamma)), 1.0f);
  return clamp(I_gamma_corrected, vec3(0.0f), vec3(1.0f));
}

// #endregion

int main(int argc, char const *argv[]) {
  clock_t t = clock(); // keeping the time of the rendering
                       // Default
                       // Final
  int width = 2048;
  int height = 1536;
  // Debug
  //  int width = 1024;
  //  int height = 768;
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
