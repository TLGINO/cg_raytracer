#ifndef RAYTRACER_MESH_H
#define RAYTRACER_MESH_H

#include "../../libs/glm/glm.hpp"
#include "hit.h"
#include "material.h"
#include "object.h"
#include "ray.h"
#include "triangle.h"
#include <fstream>
#include <iosfwd>
#include <vector>

namespace raytracer {

using std::string;
using std::vector;

class Mesh : public Object {
 public:
  Mesh(string, Material);
  Hit intersect(Ray) override;
  void load_mesh();

 private:
  vector<Triangle> triangles_;
  string file_name_;
};
}  // namespace raytracer

#endif  // RAYTRACER_MESH_H
