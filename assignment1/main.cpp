/**
@file main.cpp
@authors Jeferson Morales Mariciano, Martin Lettry
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "glm/glm.hpp"
#include "Image.h"

using glm::radians; // glm::normalize, glm::tan, glm::radians, glm::dot, glm::sqrt;
using glm::vec3;
using std::cout;
using std::endl;
using std::vector;

/**
 * Class representing a single ray.
 */
class Ray
{
public:
    vec3 origin;
    vec3 direction;

    /**
     Constructor of the ray
     @param origin Origin of the ray
     @param direction Direction of the ray
     */
    Ray(vec3 origin, vec3 direction) : origin(origin), direction(direction){};
};

class Object;

/**
 * Structure representing the event of hitting an object
 */
struct Hit
{
    bool hit = false;          // indicating whether there was or there was no intersection with an object
    vec3 normal;               // Normal vector of the intersected object at the intersection point
    vec3 intersection;         // Point of Intersection
    float distance = INFINITY; // Distance from the origin of the ray to the intersection point
    Object *object;            // A pointer to the intersected object
};

/**
 * General class for the object
 */
class Object
{
public:
    vec3 color; // Color of the object

    /**
     * A function computing an intersection, which returns the structure Hit
     */
    virtual Hit intersect(Ray ray) = 0;
};

/**
 * Implementation of the class Object for sphere shape.
 */
class Sphere : public Object
{
private:
    float radius; // Radius of the sphere
    vec3 center;  // Center of the sphere

public:
    /**
     * The constructor of the sphere
     * @param radius Radius of the sphere
     * @param center Center of the sphere
     * @param color Color of the sphere
     */
    Sphere(float radius, vec3 center, vec3 color)
        : radius(radius), center(center)
    {
        this->color = color;
    }

    /**
     * Implementation of the intersection function
     */
    Hit intersect(Ray ray)
    {
        Hit hit;
        hit.hit = false; // at start no intersection

        // ray sphere intersection. Remember to set all the fields of the hit structure:
        // hit.intersection, hit.normal, hit.distance, hit.object
        vec3 c = center - ray.origin; // vector ray - sphereCenter tracing

        // triangle ray sphere: c, ray.direction
        const float a = dot(c, ray.direction);            // length rayOrigin - sphereCenterPerpendicularNearestPoint
        const float c_norm = sqrt(dot(c, c));             // length rayOrigin - sphereCenter
        const float D = sqrt(pow(c_norm, 2) - pow(a, 2)); // length sphereCenter - sphereCenterPerpendicularNearestPoint

        // no intersection
        if (D > radius)
            return hit;

        // triangle sphere perpendicularRay: r, D
        const float b = sqrt(pow(radius, 2) - pow(D, 2));
        const float t1 = a - b, t2 = a + b;
        const float t = t1 <= 0 ? t2 : t1; // if t < 0  inside intersected object, ignore cause off camera

        // if t2 < 0 the whole intersected object is behind us
        if (t2 < 0)
            return hit;

        hit.hit = true;
        hit.object = this;
        hit.intersection = t * ray.direction + ray.origin;
        hit.distance = distance(ray.origin, hit.intersection);
        hit.normal = normalize(hit.intersection - center);
        return hit;
    }
};

vector<Object *> objects; // objects in scene

/**
 * Functions that computes a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 trace_ray(Ray ray)
{
    const vec3 COLOR_BLACK = vec3(0.0, 0.0, 0.0);

    // hit structure representing the closest intersection
    Hit closest_hit;

    // Loop over all objects to find the closest intersection
    for (int k = 0; k < objects.size(); k++)
    {
        const Hit hit = objects[k]->intersect(ray);

        if (hit.hit && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    return closest_hit.hit ? closest_hit.object->color : COLOR_BLACK;
}

/**
 * Function defining the scene
 */
void sceneDefinition()
{
    // first sphere (Exercise 1)
    objects.push_back(new Sphere(1.0, vec3(-0, -2, 8), vec3(0.6, 0.9, 0.6)));
    // additional sphere (Exercise 2)
    objects.push_back(new Sphere(1.0, vec3(1.0, -2.0, 8.0), vec3(0.6, 0.6, 0.9)));
}

int main(int argc, const char *argv[])
{
    clock_t t = clock(); // variable for keeping the time of the rendering

    // image
    const int width = 1024;
    const int height = 768;
    const float fov = 90; // field of view, horizontal opening angle, alpha

    sceneDefinition();
    Image image(width, height); // Create an image where we will store the result

    const vec3 ORIGIN(0, 0, 0);                             // camera location in front of image plane
    const float Z = 1;                                      // image plane at z
    const float size = 2 * tan(0.5 * radians(fov)) / width; // size of a pixel
    const float X = (-1) * width * size / 2;                // start of X at top left corner
    const float Y = height * size / 2;                      // start of Y at top left corner

    // Loop over pixels to form and traverse the rays through the scene
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
        {
            // ray definition for pixel (i,j), ray traversal
            const float dx = X + i * size + size / 2;
            const float dy = Y - j * size - size / 2;
            const vec3 direction(dx, dy, Z);
            const vec3 direction_norm = normalize(direction);
            const Ray ray(ORIGIN, direction_norm); // ray traversal
            image.setPixel(i, j, trace_ray(ray));
        }

    t = clock() - t;
    cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
    cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

    // Writing the final results of the rendering
    if (argc == 2)
        image.writeImage(argv[2]);
    else
        image.writeImage("./result.ppm");

    return 0;
}

/**
 * to test email about ray origin when not at point [0,0,0] and comment scene definition line at 163
    auto mysphere = Sphere(1.0, vec3(0.0, -2.0, 8.0), vec3(1.0, 0.0, 0.0));
    Ray r(glm::vec3(0, -2, 12), glm::vec3(0,0, -1));
    Hit h = mysphere.intersect(r);
    cout << h.hit << endl;
 */
