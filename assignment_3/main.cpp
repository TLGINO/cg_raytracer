/**
* @file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#include "glm/glm.hpp"
#include "Image.h"
#include "Material.h"

#define PI atan(1)*4

using glm::vec3;
using glm::normalize;
using glm::reflect;
using glm::distance;

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

class Object;

/**
 * Structure representing the event of hitting an object
 */
struct Hit {
	bool hit = false;	// state if intersection with an object
	vec3 normal;		// Normal vector of the intersected object at the intersection point
	vec3 intersection;  // Point of Intersection
	float distance = INFINITY;	// Distance from the origin of the ray to the intersection point
	Object * object = 0;	// pointer to the intersected object
	vec2 uv;	// Coordinates for computing the texture (texture coordinates)
};

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
	Material getMaterial() const {
		return material;
	}

	/**
	 * Function that set the material
	 * @param material A structure describing the material of the object
	 */
	void setMaterial(Material m) {
		this->material = m;
	}
};

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
        hit.uv.s = (asin(hit.normal.y) + M_PI/2) / M_PI;
        hit.uv.t = (atan2(hit.normal.z, hit.normal.x) + M_PI) / (2 * M_PI);
		return hit;
	}
};

class Plane : public Object {
private:
	vec3 normal;
	vec3 point;

public:
	Plane(vec3 point, vec3 normal) : point(point), normal(normal) {}
	Plane(vec3 point, vec3 normal, Material material) : point(point), normal(normal) {
		this->material = material;
	}

		float t = glm::dot(point - ray.origin, normal) / glm::dot(ray.direction, normal);

		if (t <= 0)
			return Hit();

		Hit hit;
		hit.hit = true;
		hit.intersection = ray.origin + t * ray.direction;
		hit.distance = glm::distance(ray.origin, hit.intersection);
		hit.normal = glm::normalize(normal);
		hit.object = this;
		return hit;
	}
};

/**
 * Light class
 */
class Light {
public:
	vec3 position;
	vec3 color;

    Light(vec3 position, vec3 color) : position(position), color(color) {}

	explicit Light(vec3 position) {
        Light(position, vec3(1.0f));
	}
};

// GLOBAL VARIABLES
vector<Light *> lights;     // list of lights in the scene
vec3 ambient_light(0.5f);
vector<Object *> objects;   // list of all objects in the scene

/**
 * Function for computing color of an object according to the Phong Model
 * @param point A point belonging to the object for which the color is computed
 * @param normal A normal vector the the point
 * @param uv Texture coordinates
 * @param view_direction A normalized direction from the point to the viewer/camera
 * @param material A material structure representing the material of the object
 */
vec3 PhongModel(vec3 point, vec3 normal, vec2 uv, vec3 view_direction, Material material) {
	vec3 color(0.0f);

	glm::vec3 color(0.0);
	for (int light_num = 0; light_num < lights.size(); light_num++)
	{

		glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
		glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

		float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
		float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

		glm::vec3 diffuse_color = material.diffuse;
		if (material.texture)
			diffuse_color = material.texture(uv);

		glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
		glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));

		// HERE JEFF
		/*


		 Excercise 3 - Modify the code by adding attenuation of the light due to distance from the intersection point to the light source



		 */

		color += lights[light_num]->color * (diffuse + specular);
	}
	color += ambient_light * material.ambient;

	color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
	return color;
}

/**
 * Function computing a color along the ray
 * @param ray Ray that should be traced through the scene
 * @return Color at the intersection point
 */
vec3 trace_ray(Ray ray) {
	Hit closest_hit;

	for (Object*& o : objects) {
		Hit hit = o->intersect(ray);
		if (hit.hit && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}

	return !closest_hit.hit
    ? vec3(0.0f)
    : PhongModel(closest_hit.intersection,
                 closest_hit.normal,
                 closest_hit.uv,
                 normalize(-ray.direction),
                 closest_hit.object->getMaterial());
}

/**
 * Function defining the scene
 */
void sceneDefinition()
{

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.07f, 0.09f, 0.07f);
	green_diffuse.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
	red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.07f, 0.07f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;

	Material material_rainbow;
	material_rainbow.texture = &rainbowTexture;
	material_rainbow.ambient = glm::vec3(0.0f);
	material_rainbow.specular = glm::vec3(1.0);
	material_rainbow.shininess = 10.0;

	objects.push_back(new Sphere(1.0, glm::vec3(1, -2, 8), blue_specular));
	objects.push_back(new Sphere(0.5, glm::vec3(-1, -2.5, 6), red_specular));
	objects.push_back(new Sphere(1.0, glm::vec3(3, -2, 6), green_diffuse));

	objects.push_back(new Sphere(5.0, glm::vec3(-5, 3, 20), material_rainbow));

	objects.push_back(new Plane(glm::vec3(-15, 0, 0), glm::vec3(1, 0, 0), blue_specular));
	objects.push_back(new Plane(glm::vec3(15, 0, 0), glm::vec3(-1, 0, 0), red_specular));

	objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0, -1, 0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0, 1, 0), blue_specular));

	objects.push_back(new Plane(glm::vec3(0, 0, -0.01), glm::vec3(0, 0, -1), red_specular));
	objects.push_back(new Plane(glm::vec3(0, 0, 30), glm::vec3(0, 0, 1), green_diffuse));

	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.4)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.4)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
}

/**
 * Performing tone mapping and gamma correction of intensity computed using the raytracer
 * EX3
 * @param intensity Input intensity
 * @return Tone mapped intensity in range [0,1]
 */
glm::vec3 toneMapping(glm::vec3 intensity)
{
	float gamma = 1.1f;
	float alpha = 1.0f;
	glm::vec3 tonemapped = alpha * glm::pow(intensity, glm::vec3(gamma)); // tonemapped intensity
	return glm::clamp(tonemapped, glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char *argv[]) {
	clock_t t = clock(); // keeping the time of the rendering
	int width = 1024;
	int height = 768;
	float fov = 90;	  // field of view
	sceneDefinition();
	Image image(width, height); // Create an image where we will store the result
	float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
	float X = -s * width / 2;
	float Y = s * height / 2;

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++) {
			float dx = X + i * s + s / 2;
			float dy = Y - j * s - s / 2;
			vec3 origin(0);
			vec3 direction(dx, dy, 1);
			direction = normalize(direction);
			Ray ray(origin, direction);
			image.setPixel(i, j, toneMapping(trace_ray(ray)));
		}

	t = clock() - t;
	cout << "It took " << (float)t / CLOCKS_PER_SEC << " seconds to render the image." << endl;
	cout << "I could render at " << (float)CLOCKS_PER_SEC / (float)t << " frames per second." << endl;

	// Writing the final results of the rendering
	if (argc == 2)
		image.writeImage(argv[1]);
	else
		image.writeImage("./result.ppm");

	return 0;
}
