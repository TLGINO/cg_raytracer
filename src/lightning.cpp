#include "raytracer/lightning.h"
#include <cmath>

using glm::pow;
using glm::vec2;
using glm::vec3;
using raytracer::Hit;
using raytracer::Light;
using raytracer::Material;
using raytracer::Object;
using raytracer::Ray;
using std::max;
using std::vector;

/**
 * Function for computing color of an object according to the Phong Model
 * @param point A point belonging to the object for which the color is
 * computed
 * @param normal A normal vector the the point
 * @param uv Texture coordinates
 * @param view_direction A normalized direction from the point to the
 * viewer/camera
 * @param material A material structure representing the material of the
 * object
 */
vec3 PhongModel(vec3 point, vec3 normal, vec2 uv, vec3 view_direction,
                Material material, vector<Light*> lights, vec3 ambient_light) {
  vec3 color(0.0F);

  vec3 I_diffuse = vec3(0), I_specular = vec3(0),
       // ambient illumination
      I_ambient = material.ambient * ambient_light;

  for (raytracer::Light*& l : lights) {
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
        glm::clamp(dot(normal, l_direction), 0.0f, 1.0f);  // already normalized

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
    float k = max(material.shininess, 1.0F);
    float cos_alpha = glm::clamp(dot(view_direction, r_direction), 0.0F, 1.0F);
    I_specular +=
        material.specular * pow(cos_alpha, k) * l->color * attenuation;
  }

  color = I_ambient + I_diffuse + I_specular;
  return clamp(color, vec3(0), vec3(1));
}

/**
 * Performing tone mapping and gamma correction of intensity computed using
 * the raytracer EX3 value of gamma from:
 * https://www.rtings.com/laptop/reviews/dell/precision-5560-2021#test_5194
 * @param intensity Input intensity
 * @return Tone mapped intensity in range [0,1]
 */
vec3 raytracer::lightning::toneMapping(vec3 intensity) {
  float alpha = 0.8F;
  float beta = 1.2F;
  float gamma = 2.13F;
  vec3 intensity_tone_mapped = alpha * pow(intensity, vec3(beta));
  vec3 intensity_gamma_corrected =
      min(pow(intensity_tone_mapped, vec3(1 / gamma)), 1.0F);
  return clamp(intensity_gamma_corrected, vec3(0.0F), vec3(1.0F));
}

/**
 * Function computing a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 raytracer::lightning::traceRay(vector<Object*> objects, Ray ray,
                                    vector<Light*> lights, vec3 ambient_light) {
  Hit closest_hit;

  for (Object*& o : objects) {
    Hit hit = o->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance) closest_hit = hit;
  }

  return !closest_hit.hit
             ? vec3(0.0F)
             : PhongModel(closest_hit.intersection, closest_hit.normal,
                          closest_hit.uv, normalize(-ray.direction),
                          closest_hit.object->getMaterial(), lights,
                          ambient_light);
}
