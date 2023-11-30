#include "../include/raytracer/object.h"

using raytracer::Hit;
using raytracer::Material;
using raytracer::Object;

/** Function that returns the material struct of the object*/
Material Object::getMaterial() const { return material; }

/**
 * Function that set the material
 * @param material A structure describing the material of the object
 */
[[maybe_unused]] void Object::setMaterial(Material m) { this->material = m; }
