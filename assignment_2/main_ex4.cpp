/**
 * @file main.cpp
 * @authors Jeferson Morales Mariciano, Martin Lettry
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#include "glm/glm.hpp"
#include "Image.h"
#include "Material.h"

using glm::vec3;

/**
 * Class representing a single ray.
 */
class Ray{
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

class Object;

/**
 * Structure representing the even of hitting an object
 */
struct Hit {
    bool hit = false; // state intersection with an object
    vec3 normal; // Normal vector of the intersected object at the intersection point
    vec3 intersection; // Point of Intersection
    float distance = INFINITY; // Distance from the origin of the ray to the intersection point
    Object * object; // pointer to the intersected object
};

/**
 * General class for the object
 */
class Object {
public:
    vec3 color = vec3(0);
    Material material; // Structure describing the material of the object

    /** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

    /** Function that returns the material struct of the object */
    Material getMaterial() const {
        return material;
    }

    /**
     * Function that set the material
     * @param material structure describing the material of the object
    */
    void setMaterial(Material m){
        this->material = m;
    }
};

/**
 * Implementation of the class Object for sphere shape.
 */
class Sphere : public Object{
private:
    float radius;
    vec3 center;

public:
    /**
     * The constructor of the sphere
     * @param radius Radius of the sphere
     * @param center Center of the sphere
     * @param color Color of the sphere
     */
    Sphere(float radius, vec3 center, vec3 color) : radius(radius), center(center) {
        this->color = color;
    }

    Sphere(float radius, vec3 center, Material material) : radius(radius), center(center) {
        this->material = material;
    }

    /** Implementation of the intersection function */
    Hit intersect(Ray ray) override {
        vec3 c = center - ray.origin;

        float cdotc = dot(c,c);
        float cdotd = dot(c, ray.direction);

        Hit hit;

        float D = 0;
        if (cdotc > cdotd*cdotd) {
            D =  sqrt(cdotc - cdotd*cdotd);
        }

        if (D > radius) {
            hit.hit = false;
            return hit;
        }

        hit.hit = true;
        float t1 = cdotd - sqrt(radius*radius - D*D);
        float t2 = cdotd + sqrt(radius*radius - D*D);

        float t = t1;
        if(t < 0) t = t2;
        if(t < 0) {
            hit.hit = false;
            return hit;
        }

        hit.intersection = ray.origin + t * ray.direction;
        hit.normal = normalize(hit.intersection - center);
        hit.distance = distance(ray.origin, hit.intersection);
        hit.object = this;
        return hit;
    }
};

/**
 * Light class
 */
class Light{
public:
    vec3 position;
    vec3 color;

    explicit Light(vec3 position): position(position), color(vec3(1)) {}
    Light(vec3 position, vec3 color): position(position), color(color) {}
};

// GLOBAL VARIABLES
vector<Light *> lights; // list of lights in the scene
vec3 ambient_light(1);
vector<Object *> objects; // list of all objects in the scene

/**
 * Function for computing color of an object according to the Phong Model
 * @param point A point belonging to the object for which the color is computed
 * @param normal normal vector from the point
 * @param view_direction A normalized direction from the point to the viewer/camera
 * @param material A material structure representing the material of the object
*/
vec3 PhongModel(vec3 point, vec3 normal, vec3 view_direction, Material material) {
    vec3 color(0);
    auto I_diffuse = vec3(0), I_specular = vec3(0);
    for (auto & light : lights) {
        // diffuse reflection: from diffuse reflection point, cos of angle of surface normal n and direction from point to light source
        auto l = normalize(light->position - point);
        auto phi = max(dot(normal,l), 0.0f);  // normal already normalize
        I_diffuse += material.diffuse * phi * light->color;

        // specular highlight
        auto r = reflect(-l, normal);
        auto k = max(material.shininess, 1.0f);
        auto cos_alpha = max(dot(view_direction, r), 0.0f);
        I_specular += material.specular * pow(cos_alpha, k) * light->color;
    }
    // ambient illumination
    auto I_ambient = material.ambient * ambient_light;

    color = I_ambient + I_diffuse + I_specular;
    // The final color has to be CLAMPED so the values do not go beyond 0 and 1.
    return clamp(color, vec3(0), vec3(1));
}

/**
 * Functions that computes a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 trace_ray(Ray ray){
    Hit closest_hit;

    for (auto & object : objects) {
        Hit hit = object->intersect(ray);
        if (hit.hit && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    return !closest_hit.hit
           ? vec3(0)  // black
           : PhongModel(closest_hit.intersection,
                        closest_hit.normal,
                        normalize(-ray.direction),
                        closest_hit.object->getMaterial());
}

/**
 * Function defining the scene
 */
void sceneDefinition () {
    auto sphere_red = new Sphere(0.5f, {-1, -2.5, 6}, {1, 0, 0});
    auto sphere_blue = new Sphere(1,{1, -2, 8},{0, 0, 1});
    auto sphere_green = new Sphere(1, {3, -2, 6}, {0, 1, 0});

    sphere_red->setMaterial(Material(
            {.ambient = {0.01, 0.03, 0.03}, .diffuse = {1, 0.3, 0.3}, .specular = vec3(0.5), .shininess = 10}));
    sphere_blue->setMaterial(Material(
            {.ambient = {0.07, 0.07, 0.1}, .diffuse = {0.7, 0.7, 1}, .specular = vec3(0.6), .shininess = 100.0}));
    sphere_green->setMaterial(Material(
            {.ambient = {0.07, 0.09, 0.07}, .diffuse = {0.7, 0.9, 0.7}, .specular = vec3(0), .shininess = 0}));

    objects.push_back(sphere_red);
    objects.push_back(sphere_blue);
    objects.push_back(sphere_green);

    lights.push_back(new Light({0, 26, 5}, vec3(0.4)));
    lights.push_back(new Light({0, 1, 12}, vec3(0.4)));
    lights.push_back(new Light({0, 5, 1}, vec3(0.4)));
}

int main(int argc, const char *argv[]) {
    clock_t t = clock(); // keeping the time of the rendering
    const int width = 1024;
    const int height = 768;
    const float fov = 90; // field of view
    sceneDefinition();

    if (argc == 2)
    {
        try
        {
            float angle = 2 * M_PI / 60.0 * stoi(argv[1]);
            printf("Rendering frame with a light at an angle of %f \n", angle);

            float light_x = lights[1]->position.x * cos(angle) - lights[1]->position.z * sin(angle);
            float light_z = lights[1]->position.x * sin(angle) + lights[1]->position.z * cos(angle);
            lights[1]->position = glm::vec3(light_x, lights[1]->position.y, light_z);
        }
        catch (...)
        {
            printf("Invalid input given for angle, expected a float, got %s", argv[2]);
        }
    }

    Image image(width, height); // Create an image where we will store the result
    float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
    float X = -s * width / 2;
    float Y = s * height / 2;
    const float Z = 1;

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {

            float dx = X + i * s + s / 2;
            float dy = Y - j * s - s / 2;

            glm::vec3 origin(0);
            glm::vec3 direction(dx, dy, Z);
            direction = normalize(direction);
            Ray ray(origin, direction);
            image.setPixel(i, j, trace_ray(ray));
        }

    t = clock() - t;
    cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
    cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

    // Writing the final results of the rendering
    if (argc == 2) {
        std::string f_name = std::string(argv[1]) + ".ppm";
        image.writeImage(f_name.c_str());
    } else
        image.writeImage("./result.ppm");

    return 0;
}

/*
 * SHELL SCRIPT TO GENERATE VIDEO
#!/bin/bash

g++ main.cpp -o a.out

for i in $(seq 1 180); do
  ./a.out $i
done

input_pattern="%01d.ppm" # Change this pattern to match your file naming

output_file="output.mp4"
ffmpeg -framerate 60 -i "$input_pattern" -c:v vp9  -r 60 "$output_file" -y

rm [0-9]*.ppm
*/
