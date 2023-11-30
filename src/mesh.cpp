#include "../include/raytracer/mesh.h"
#include <iostream>
#include <sstream>
#include <string>

using raytracer::Hit;
using raytracer::Mesh;
using std::endl;
using std::string;
using std::vector;

Mesh::Mesh(std::string file_name, Material material) {
  this->material = material;
  this->file_name_ = file_name;
  this->load_mesh();
}

void Mesh::load_mesh() {
  std::ifstream objFile(this->file_name_);
  std::string line;
  vector<vec3> v_positions;
  vector<vec3> vn_positions;

  while (getline(objFile, line)) {
    std::istringstream stream(line);
    std::string token;
    stream >> token;

    if (token == "v") {
      float x, y, z;
      stream >> x >> y >> z;
      v_positions.push_back(glm::vec3(x, y, z));
    }
    if (token == "vn") {
      float x, y, z;
      stream >> x >> y >> z;
      vn_positions.push_back(glm::vec3(x, y, z));
    }
    if (token == "f") {
      string i, j, k;
      stream >> i >> j >> k;

      if (i.find("//") != std::string::npos) {
        int i_v = stoi(i.substr(0, i.find("//")));
        int i_vn = stoi(i.substr(i.find("//") + 2, i.size()));

        int j_v = stoi(j.substr(0, j.find("//")));
        int j_vn = stoi(j.substr(j.find("//") + 2, j.size()));

        int k_v = stoi(k.substr(0, k.find("//")));
        int k_vn = stoi(k.substr(k.find("//") + 2, k.size()));

        Triangle t = Triangle(v_positions[i_v - 1], v_positions[j_v - 1],
                              v_positions[k_v - 1], vn_positions[i_vn - 1],
                              vn_positions[j_vn - 1], vn_positions[k_vn - 1],
                              this->material);
        this->triangles_.push_back(t);
      } else {
        Triangle t =
            Triangle(v_positions[stoi(i) - 1], v_positions[stoi(j) - 1],
                     v_positions[stoi(k) - 1], this->material);
        this->triangles_.push_back(t);
      }
    }
  }
  objFile.close();
}

Hit Mesh::intersect(Ray ray) {
  // Does this logic make sense?
  for (Triangle t : this->triangles_) {
    Hit h = t.intersect(ray);
    if (h.hit) {
      std::cout << "HIT MESH" << endl;
      return h;
    }
  }
  return Hit();
}
