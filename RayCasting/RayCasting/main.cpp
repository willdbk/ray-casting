

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <time.h>


#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

#include <vector>

#include "vec3.h"
#include "vec2.h"
#include "vec4.h"
#include "mat4x4.h"

const unsigned int windowWidth = 512, windowHeight = 512;

int majorVersion = 3, minorVersion = 0;

void getErrorInfo(unsigned int handle)
{
    int logLen;
    glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &logLen);
    if (logLen > 0)
    {
        char * log = new char[logLen];
        int written;
        glGetShaderInfoLog(handle, logLen, &written, log);
        printf("Shader log:\n%s", log);
        delete log;
    }
}

void checkShader(unsigned int shader, char * message)
{
    int OK;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &OK);
    if (!OK)
    {
        printf("%s!\n", message);
        getErrorInfo(shader);
    }
}

void checkLinking(unsigned int program)
{
    int OK;
    glGetProgramiv(program, GL_LINK_STATUS, &OK);
    if (!OK)
    {
        printf("Failed to link shader program!\n");
        getErrorInfo(program);
    }
}

class Shader
{
protected:
    unsigned int shaderProgram;

public:
    Shader()
    {
        const char *vertexSource = "\n\
        #version 150 \n\
        precision highp float; \n\
        \n\
        in vec2 vertexPosition;	\n\
        in vec2 vertexTexCoord; \n\
        out vec2 texCoord; \n\
        \n\
        void main() \n\
        { \n\
        texCoord = vertexTexCoord; \n\
        gl_Position = vec4(vertexPosition.x, vertexPosition.y, 0, 1); \n\
        } \n\
        ";

        const char *fragmentSource = "\n\
        #version 150 \n\
        precision highp float; \n\
        \n\
        uniform sampler2D samplerUnit; \n\
        in vec2 texCoord;  \n\
        out vec4 fragmentColor; \n\
        \n\
        void main() { \n\
        fragmentColor = texture(samplerUnit, texCoord);  \n\
        } \n\
        ";

        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        if (!vertexShader) { printf("Error in vertex shader creation\n"); exit(1); }

        glShaderSource(vertexShader, 1, &vertexSource, NULL);
        glCompileShader(vertexShader);
        checkShader(vertexShader, "Vertex shader error");

        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fragmentShader) { printf("Error in fragment shader creation\n"); exit(1); }

        glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
        glCompileShader(fragmentShader);
        checkShader(fragmentShader, "Fragment shader error");

        shaderProgram = glCreateProgram();
        if (!shaderProgram) { printf("Error in shader program creation\n"); exit(1); }

        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);

        glBindAttribLocation(shaderProgram, 0, "vertexPosition");
        glBindAttribLocation(shaderProgram, 1, "vertexTexCoord");

        glBindFragDataLocation(shaderProgram, 0, "fragmentColor");

        glLinkProgram(shaderProgram);
        checkLinking(shaderProgram);
    }

    ~Shader()
    {
        if(shaderProgram) glDeleteProgram(shaderProgram);
    }

    void Run()
    {
        if(shaderProgram) glUseProgram(shaderProgram);
    }

    void UploadSamplerID()
    {
        int samplerUnit = 0;
        int location = glGetUniformLocation(shaderProgram, "samplerUnit");
        glUniform1i(location, samplerUnit);
        glActiveTexture(GL_TEXTURE0 + samplerUnit);
    }
};

Shader *shader = 0;

// Simple material class, with object color, and headlight shading.

//class Material
//{
//    vec3 frontColor, backColor;
//public:
//    Material(vec3 frontColor = vec3(1,1,1), vec3 backColor = vec3(1,1,1)): frontColor(frontColor), backColor(backColor) {
//    }
//
//    virtual vec3 getColor(vec3 position, vec3 normal, vec3 viewDir)
//    {
//        double w = normal.dot(viewDir);
//        if(w < 0) {
//            return backColor*(-w);
//        }
//        return frontColor*w;
//    }
//};

class Material
{
    vec3 color;
    vec3 ka, kd, ks;
    float shininess;

public:
    bool pinked;

    Material(vec3 color = vec3(1,1,1), bool pinked = false, vec3 ka = vec3(0.2,0.2,0.2), vec3 kd = vec3(0.6, 0.6, 0.6), vec3 ks = vec3(0.3, 0.3, 0.3), float shininess = 50.0): color(color), pinked(pinked), ka(ka), kd(kd), ks(ks), shininess(shininess) {}

    virtual vec3 getColor(vec3 position, vec3 normal, vec3 viewDir) {
        return color;
    }

    virtual vec3 shade(vec3 position, vec3 normal, vec3 viewDir, vec3 lightDir, vec3 powerDensity)
    {
        vec3 N = normal.normalize();
        vec3 V = viewDir.normalize();
        vec3 L = lightDir.normalize();
        vec3 H = (V + L).normalize();
        vec3 La = vec3(0,0,0);
        vec3 raw_color =
            La * ka +
            powerDensity * kd * getColor(position, normal, viewDir) * fmax(0.0, L.dot(N)) +
            powerDensity * ks * pow(fmax(0.0, H.dot(N)), shininess);
        return raw_color;
    }
};



float snoise(vec3 r) {
    unsigned int x = 0x0625DF73;
    unsigned int y = 0xD1B84B45;
    unsigned int z = 0x152AD8D0;
    float f = 0;
    for(int i=0; i<32; i++) {
        vec3 s(	x/(float)0xffffffff,
               y/(float)0xffffffff,
               z/(float)0xffffffff);
        f += sin(s.dot(r));
        x = x << 1 | x >> 31;
        y = y << 1 | y >> 31;
        z = z << 1 | z >> 31;
    }
    return f / 64.0 + 0.5;
};


class Wood : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Wood(): Material(vec3(1, 1, 1))
    {
        scale = 16;
        turbulence = 500;
        period = 8;
        sharpness = 10;
    }
    virtual vec3 getColor(vec3 position, vec3 normal, vec3 viewDir)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
        return (vec3(1, 0.3, 0) * w + vec3(0.35, 0.1, 0.05) * (1-w));
    }
};

class Marble : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Marble():
    Material(vec3(1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 32;
        sharpness = 1;
    }
    virtual vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        //return normal;
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        return (vec3(0, 0, 1) * w + vec3(1, 1, 1) * (1-w)) * normal.dot(viewDir);
    }
};


class Fire : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Fire(bool pinked = false): Material(vec3(1, 1, 1), pinked)
    {
        scale = 16;
        turbulence = 500;
        period = 8;
        sharpness = 10;
    }
    virtual vec3 getColor(vec3 position, vec3 normal, vec3 viewDir)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
        if(w < 0.3) { return vec3(244, 132, 22).colorize();}
        else if(w < 0.4) { return vec3(243, 203, 33).colorize();}
        else if(w < 0.5) { return vec3(246, 226, 33).colorize();}
        else if(w < 0.6) { return vec3(248, 245, 88).colorize();}
        return vec3(224, 56, 10).colorize();

    }
    virtual vec3 shade(vec3 position, vec3 normal, vec3 viewDir, vec3 lightDir, vec3 powerDensity)
    {
        return getColor(position, normal, viewDir);
    }
};


class Aurora : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Aurora(): Material(vec3(1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 3;
        sharpness = 1;
    }
    virtual vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        //return normal;
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        w = w - int(w)-0.5;

        float a = (0.8*position.x + 1.2*position.y) + 2*w;
        if(sin((position.x-position.y)*period+ snoise(position * scale))> 0.4) {
            a += sin((position.x-position.y)*period+ snoise(position * scale))-0.4;
        }
        if(sin((position.x-position.y)*period+ snoise(position * scale)) < -0.4) {
            a += sin((position.x-position.y)*period+ snoise(position * scale))+0.4;
        }
        if(a <= -2) { return vec3(20,232,30).colorize(); }
        if(a <= 0) { return vec3(0,234,141).colorize(); }
        if(a <= 2) { return vec3(1,126,213).colorize(); }
        if(a <= 4) { return vec3(181,61,255).colorize(); }
        else { return vec3(141,0,196).colorize(); }


        return (vec3(0, 0, 1) * w + vec3(1, 1, 1) * (1-w)) * normal.dot(viewDir);


    }

};


class Sky
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Sky()
    {
        scale = 10;
        turbulence = 200;
        period = 8;
        sharpness = 10;
    }
    vec3 getColor(vec3 position)
    {
        float w = pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);

        return vec3(135/255.0,206/255.0,250/255.0)*w+vec3(1,1,1)*(1-w);
    }
};


// Camera class.

class Camera
{
    vec3 eye;		// World space camera position.
    vec3 lookAt;	// Center of window in world space.
    vec3 right;		// Vector from window center to window right-mid (in world space).
    vec3 up;		// Vector from window center to window top-mid (in world space).

public:
    Camera()
    {
        eye = vec3(0, 0, 2);
        lookAt = vec3(0, 0, 1);
        right = vec3(1, 0, 0);
        up = vec3(0, 1, 0);
    }
    vec3 getEye()
    {
        return eye;
    }

    // Compute ray through pixel at normalized device coordinates.

    vec3 rayDirFromNdc(float x, float y) {
        return (lookAt - eye
                + right * x
                + up    * y
                ).normalize();
    }
};

// Light Hierarchy
class LightSource
{
public:
    virtual vec3 getPowerDensityAt ( vec3 x )=0;
    virtual vec3 getLightDirAt     ( vec3 x )=0;
    virtual float  getDistanceFrom ( vec3 x )=0;
};

class DirectionalLight : public LightSource
{
public:
    vec3 power;
    vec3 dir;

    DirectionalLight(vec3 dir, vec3 power): dir(dir), power(power) {}

    vec3 getPowerDensityAt ( vec3 x ) {
        return power;
    }
    vec3 getLightDirAt     ( vec3 x ) {
        return dir;
    }
    float getDistanceFrom ( vec3 x ) {
        return FLT_MAX;
    }
};


class PointLight : public LightSource
{
    vec3 power;
    vec3 pos;
public:
    PointLight(vec3 pos, vec3 power): pos(pos), power(power) {}

    vec3 getPowerDensityAt ( vec3 x ) {
        return power*(1/pow(getDistanceFrom(x),2));
    }
    vec3 getLightDirAt     ( vec3 x ) {
        return (pos-x).normalize();
    }
    float  getDistanceFrom ( vec3 x ) {
        float dist = (x-pos).norm();
        return dist;
    }
};



// Ray structure.

class Ray
{
public:
    vec3 origin;
    vec3 dir;
    bool isShadowRay;
    Ray(vec3 o, vec3 d, bool isShadowRay = false): isShadowRay(isShadowRay)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.

class Hit
{
public:
    Hit()
    {
        t = -1;
    }
    float t;				// Ray paramter at intersection. Negative means no valid intersection.
    vec3 position;			// Intersection coordinates.
    vec3 normal;			// Surface normal at intersection.
    Material* material;		// Material of intersected surface.
};

// Abstract base class.

class Intersectable
{
protected:
    Material* material;
    bool shadowCaster;
public:
    Intersectable(Material* material, bool _shadowCaster = true):material(material) {
        shadowCaster = _shadowCaster;
    }
    virtual Hit intersect(const Ray& ray)=0;
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.

class QuadraticRoots
{
public:
    float t1;
    float t2;

    // Solves the quadratic a*t*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and sets members t1 and t2 to store the roots.

    QuadraticRoots(float a, float b, float c)
    {
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
        {
            t1 = -1;
            t2 = -1;
            return;
        }
        float sqrt_discr = sqrt( discr );
        t1 = (-b + sqrt_discr)/2.0/a;
        t2 = (-b - sqrt_discr)/2.0/a;
    }

    // Returns the lesser of the positive solutions, or a negative value if there was no positive solution.

    float getLesserPositive()
    {
        return (0 < t1 && (t2 < 0 || t1 < t2)) ? t1 : t2;
    }
};


// CLASS PLANE COULD COME HERE
class Plane : public Intersectable
{
    vec3 n;
    vec3 r0;
public:
    Plane(vec3 n, vec3 r0, Material* material): Intersectable(material), n(n), r0(r0) { }


    vec3 getNormalAt(vec3 r)
    {
        return n;
    }

    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.

        Hit hit;
        hit.t = -1;
        if(shadowCaster || !ray.isShadowRay) {
            float denom = ray.dir.dot(n);
            float t = -1;
            if(denom != 0) {
                t = (r0-ray.origin).dot(n)/denom;
            }

            hit.t = t;
            hit.material = material;
            hit.position = ray.origin + ray.dir * t;
            hit.normal = getNormalAt(hit.position);
        }
        return hit;
    }
};


// CLASS QUADRIC COULD COME HERE
class Quadric : public Intersectable
{
    mat4x4 coeffs;
public:
    Quadric(Material* material = 0, bool shadowCaster = true, mat4x4 coeffs = mat4x4()): Intersectable(material, shadowCaster), coeffs(coeffs) { }

    Quadric* transform(mat4x4 t) {
        coeffs = t.invert()*coeffs*t.invert().transpose();
        return this;
    }

    Quadric* ellipsoid() {
        coeffs._33 = -1;
        return this;
    }

    Quadric* cylinder() {
        coeffs._11 = 0;
        coeffs._33 = -1;
        return this;
    }

    Quadric* inverted_cylinder() {
        coeffs._00 = -1;
        coeffs._11 = 0;
        coeffs._22 = -1;
        return this;
    }

    Quadric* cone() {
        coeffs._11 = -1;
        coeffs._33 = 0;
        return this;
    }

    Quadric* paraboloid() {
        coeffs._11 = 0;
        coeffs._13 = -1;
        coeffs._33 = 0;
        return this;
    }
    Quadric* hyperboloid() {
        coeffs._11 = -1;
        coeffs._33 = -1;
        return this;
    }

    Quadric* parallelPlanes() {
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -1;
        return this;
    }

    Quadric* plane() {
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = 0;
        return this;
    }


    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        vec4 origin = vec4(ray.origin);
        vec4 dir = vec4(ray.dir.x,ray.dir.y,ray.dir.z,0);
        float a = dir.dot(coeffs*dir);
        float b = dir.dot(coeffs*origin) + origin.dot(coeffs*dir);
        float c = origin.dot(coeffs*origin);
        return QuadraticRoots(a, b, c);
    }

    vec3 getNormalAt(vec3 r)
    {
        vec4 normal = coeffs*vec4(r) + vec4(r)*coeffs;
        return vec3(normal.x,normal.y,normal.z).normalize();
    }

    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        Hit hit;
        hit.t = -1;
        if(shadowCaster || !ray.isShadowRay) {
            float t = solveQuadratic(ray).getLesserPositive();
            hit.t = t;
            hit.material = material;
            hit.position = ray.origin + ray.dir * t;
            hit.normal = getNormalAt(hit.position);
        }
        return hit;
    }

    bool contains(vec3 r) {
        vec4 rhomo(r);
        float result = (rhomo*coeffs).dot(rhomo);
        return result >=0;
    }
};



// CLASS CLIPPEDQUADRIC COULD COME HERE
class ClippedQuadric: public Intersectable {
    Quadric* shape;
    std::vector<Quadric*> clippers;
public:
    ClippedQuadric(Material* material = 0, bool shadowCaster = true, Quadric* shape = 0, Quadric* clipper = 0): Intersectable(material, shadowCaster), shape(shape) {
        if(clipper != 0) {
            clippers.push_back(clipper);
        }
    }

    ClippedQuadric* cylinder(float height) {
        shape = (new Quadric())->cylinder();
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::scaling(vec3(0.5*height,0.5*height,0.5*height))*mat4x4::translation(vec3(0,0.5*height,0))));
        return this;
    }

    ClippedQuadric* cone(float height) {
        shape = (new Quadric())->cone();
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::scaling(vec3(0.5*height,0.5*height,0.5*height))*mat4x4::translation(vec3(0,0.5*height,0))));
        return this;
    }

    ClippedQuadric* paraboloid(float height) {
        shape = (new Quadric())->paraboloid();
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::scaling(vec3(0.5*height,0.5*height,0.5*height))*mat4x4::translation(vec3(0,0.5*height,0))));
        return this;
    }

    ClippedQuadric* hyperboloid(float height) {
        shape = (new Quadric())->hyperboloid();
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::scaling(vec3(0.5*height,0.5*height,0.5*height))*mat4x4::translation(vec3(0,0.5*height,0))));
        return this;
    }

    ClippedQuadric* cube() {
        shape = (new Quadric())->parallelPlanes();
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::rotation(vec3(1.0, 0, 0), M_PI/2)));
        clippers.push_back((new Quadric())->parallelPlanes()->transform(mat4x4::rotation(vec3(0, 0, 1.0), M_PI/2)));
        return this;
    }

    ClippedQuadric* disc(float radius) {
        shape = (new Quadric())->plane();
        clippers.push_back((new Quadric())->cylinder()->transform(mat4x4::scaling(vec3(radius,1.0,radius))));
        return this;
    }

    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        Hit hit;
        hit.t = -1;
        if((shadowCaster || !ray.isShadowRay)) {
            float t1 = shape->solveQuadratic(ray).t1;
            float t2 = shape->solveQuadratic(ray).t2;
            for(Quadric* clipper : clippers) {
                if(clipper->contains(ray.origin + ray.dir * t1)) { t1 = -1; }
                if(clipper->contains(ray.origin + ray.dir * t2)) { t2 = -1; }
            }

            float t = (0 < t1 && (t2 < 0 || t1 < t2)) ? t1 : t2;

            vec3 pos = ray.origin + ray.dir * t;

            float noise = snoise(pos*100);

            if(!material->pinked || noise < 0.55) {
                hit.t = t;
                hit.material = material;
                hit.position = pos;
                hit.normal = shape->getNormalAt(hit.position);
            }
        }
        return hit;
    }

    ClippedQuadric* transform(mat4x4 t) {
        shape->transform(t);
        for(Quadric* clipper : clippers) {
            clipper->transform(t);
        }
        return this;
    }
};


class Cluster {

public:
    std::vector<Quadric*> quadrics;
    std::vector<ClippedQuadric*> clippedQuadrics;


    void transform(mat4x4 t) {
        for(Quadric* quadric: quadrics) {
            quadric->transform(t);
        }
        for(ClippedQuadric* clippedQuadric: clippedQuadrics) {
            clippedQuadric->transform(t);
        }
    }

};


class Scene
{
    Camera camera;
    std::vector<Intersectable*> objects;
    std::vector<Cluster*> clusters;
    std::vector<Material*> materials;
    std::vector<LightSource*> lightSources;
    Sky background;
public:
    Scene(Sky background = Sky()): background(background)
    {
        srand(time(NULL));

        lightSources.push_back(new DirectionalLight(vec3(0,10.0,10.0), vec3(0.5, 0.5, 0.5)));
        lightSources.push_back(new PointLight(vec3(-0.8,0.5,2.0), vec3(2.0, 2.0, 2.0)));

        Cluster* box = new Cluster();
        materials.push_back(new Marble());
        box->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->cube());
        box->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->cube()->transform(mat4x4::rotation(vec3(1,0,0), M_PI/2)));
        box->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->cube()->transform(mat4x4::rotation(vec3(0,0,1), M_PI/2)));

        box->transform(mat4x4::scaling(vec3(0.25,0.25,0.25))*mat4x4::rotation(vec3(0,1,0), -M_PI/8)*mat4x4::translation(vec3(-0.5,-0.75,0)));

        clusters.push_back(box);

        Cluster* candle = new Cluster();
        materials.push_back(new Material(vec3(.9373,.902,.8275)));
        candle->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->cylinder(1.0)->transform(mat4x4::scaling(vec3(0.25,1.0,0.25))));
        candle->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->disc(0.25)->transform(mat4x4::translation(vec3(0,0.0,0))));
        candle->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->disc(0.25)->transform(mat4x4::translation(vec3(0,1.0,0))));

        materials.push_back(new Fire(true));
        candle->quadrics.push_back((new Quadric(materials[materials.size()-1], false))->ellipsoid()->transform(mat4x4::scaling(vec3(0.05,.2,0.05))*mat4x4::translation(vec3(0,1.0,0))));
        candle->transform(mat4x4::scaling(vec3(0.5,0.5,0.5))*mat4x4::translation(vec3(0.5,-0.25,0)));
        lightSources.push_back(new PointLight(vec3(0.5,0.7,0.0), vec3(1, 1, 1)));


        clusters.push_back(candle);


        Cluster* bell = new Cluster();
        materials.push_back(new Material(vec3(1.0, 215.0/255, 0)));
        bell->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->paraboloid(1.0)->transform(mat4x4::scaling(vec3(0.5, 1.0, 0.5))));
        bell->clippedQuadrics.push_back((new ClippedQuadric(materials[materials.size()-1]))->hyperboloid(1.0)->transform(mat4x4::scaling(vec3(0.5, 1.0, 0.5))*mat4x4::translation(vec3(0, 1.0, 0))));
        bell->quadrics.push_back((new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.25, 0.5, 0.25))*mat4x4::translation(vec3(0, 2.0, 0))));

        bell->transform(mat4x4::rotation(vec3(1,0,0), M_PI)*mat4x4::scaling(vec3(0.25,0.25,0.25))*mat4x4::translation(vec3(0.5, 1.0, 0)));
        clusters.push_back(bell);



        Cluster* oranges = new Cluster();
        materials.push_back(new Material(vec3(1, 140.0/255, 0)));
        for(int i = 1; i < 5; i++) {
            for(int x = 0; x < i; x++) {
                for(int z = 0; z < i; z++) {
                    Quadric* orange = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(.1, .1, .1))*mat4x4::translation(vec3(x*(-0.2)-(4-i)*0.1,i*(-0.15),z*(0.2)+(4-i)*0.1)));
                    oranges->quadrics.push_back(orange);
                }
            }
        }
        oranges->transform(mat4x4::translation(vec3(2.8,-0.3,-2)));
        clusters.push_back(oranges);


        Cluster* snowman = new Cluster();
        materials.push_back(new Material(vec3(1, 250.0/255, 250.0/255)));

        ClippedQuadric* base = (new ClippedQuadric(materials[materials.size()-1], true, (new Quadric())->ellipsoid(), (new Quadric())->inverted_cylinder()->transform(mat4x4::scaling(vec3(0.5, 0.5, 0.5))*mat4x4::rotation(vec3(1,0,0), 90*M_PI/180))))->transform(mat4x4::scaling(vec3(0.6, 0.5, 0.6)));
        snowman->clippedQuadrics.push_back(base);
        Quadric* body = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.48, 0.4, 0.48))*mat4x4::translation(vec3(0,0.72,0)));
        snowman->quadrics.push_back(body);
        body = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.36, 0.3, 0.36))*mat4x4::translation(vec3(0,1.28,0)));
        snowman->quadrics.push_back(body);

        materials.push_back(new Material(vec3(6.0/255, 6.0/255, 7.0/255)));
        for(int i = -1; i < 2; i++) {
            Quadric* button = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.03, 0.03, 0.03))*mat4x4::translation(vec3(0,0.72+0.15*i,0.49-0.03*abs(i))));
            snowman->quadrics.push_back(button);
        }
        Quadric* eye = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.03, 0.03, 0.03))*mat4x4::translation(vec3(0.1,1.3,0.39)));
        snowman->quadrics.push_back(eye);
        eye = (new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::scaling(vec3(0.03, 0.03, 0.03))*mat4x4::translation(vec3(-0.1,1.3,0.39)));
        snowman->quadrics.push_back(eye);
        materials.push_back(new Material(vec3(0.929, 0.569, 0.129)));
        ClippedQuadric* carrot = (new ClippedQuadric(materials[materials.size()-1]))->cone(0.2)->transform(mat4x4::scaling(vec3(0.2, 1.0, 0.2))*mat4x4::rotation(vec3(1,0,0), -90.0*M_PI/180)*mat4x4::translation(vec3(0,1.26,0.55)));
        snowman->clippedQuadrics.push_back(carrot);
        snowman->transform(mat4x4::scaling(vec3(1,1,1))*mat4x4::rotation(vec3(0, 1, 0), 40.0*M_PI/180)*mat4x4::translation(vec3(-2,-0.5,-2)));
        clusters.push_back(snowman);


        Cluster* tree = new Cluster();
        materials.push_back(new Material(vec3(34/255.0,139/255.0,34/255.0), true));
        ClippedQuadric* leaves = (new ClippedQuadric(materials[materials.size()-1]))->cone(1.0)->transform(mat4x4::scaling(vec3(0.2, 0.6, 0.2))*mat4x4::rotation(vec3(1,0,0), M_PI)*mat4x4::translation(vec3(0,0.85,0)));
        tree->clippedQuadrics.push_back(leaves);
        leaves = (new ClippedQuadric(materials[materials.size()-1]))->cone(1.0)->transform(mat4x4::scaling(vec3(0.18, 0.5, 0.18))*mat4x4::rotation(vec3(1,0,0), M_PI)*mat4x4::translation(vec3(0,1.0,0)));
        tree->clippedQuadrics.push_back(leaves);
        leaves = (new ClippedQuadric(materials[materials.size()-1]))->cone(1.0)->transform(mat4x4::scaling(vec3(0.15, 0.4, 0.15))*mat4x4::rotation(vec3(1,0,0), M_PI)*mat4x4::translation(vec3(0,1.15,0)));
        tree->clippedQuadrics.push_back(leaves);
        materials.push_back(new Material(vec3(.396, .263, .129)));
        ClippedQuadric* trunk = (new ClippedQuadric(materials[materials.size()-1]))->cylinder(1.0)->transform(mat4x4::scaling(vec3(0.05, 1.0, 0.05)));
        tree->clippedQuadrics.push_back(trunk);

        tree->transform(mat4x4::scaling(vec3(2,2,2))*mat4x4::translation(vec3(-0.5,-1.0,-1)));
        clusters.push_back(tree);


        // floor
        materials.push_back(new Wood());
        objects.push_back(new Plane(vec3(0,1,0), vec3(0,-1,0), materials[materials.size()-1]));


        // side wall
        materials.push_back(new Material(vec3(1.0,1.0,224.0/255)));
        objects.push_back(new Plane(vec3(1,0,0), vec3(-3,0,0), materials[materials.size()-1]));

        // side wall
        objects.push_back(new Plane(vec3(-1,0,0), vec3(3,0,0), materials[materials.size()-1]));

        // back wall
        materials.push_back(new Aurora());
        objects.push_back(new Plane(vec3(0,0,1), vec3(0,0,-3), materials[materials.size()-1]));


        for(Cluster* cluster: clusters) {
            for(Quadric* quadric: cluster->quadrics) {
                objects.push_back(quadric);
            }

            for(ClippedQuadric* clippedQuadric: cluster->clippedQuadrics) {
                objects.push_back(clippedQuadric);
            }
        }

    }
    ~Scene()
    {
        // UNCOMMENT THESE WHEN APPROPRIATE
        for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
        	delete *iMaterial;
        for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
        	delete *iObject;
        for (std::vector<Cluster*>::iterator iCluster = clusters.begin(); iCluster != clusters.end(); ++iCluster)
            delete *iCluster;
    }

public:
    Camera& getCamera()
    {
        return camera;
    }

    Hit firstIntersect(const Ray& ray) {
        Hit bestHit;
        bestHit.t = FLT_MAX;
        for(Intersectable* obj : objects)
        {
            Hit hit = obj->intersect(ray);
            if(hit.t > 0 && hit.t < bestHit.t)
                bestHit = hit;
        }
        return bestHit;
    }



    vec3 trace(const Ray& ray)
    {
        Hit hit = firstIntersect(ray);
        if(hit.t < FLT_MAX) {
            vec3 color_sum = vec3(0.1,0.1,0.1);
            for(LightSource* source: lightSources) {
                Ray toLight = Ray(hit.position, source->getLightDirAt(hit.position));
                toLight.isShadowRay = true;
                Hit occluder = firstIntersect(toLight);
                vec3 power = source->getPowerDensityAt(hit.position);
                if(occluder.t > 0.00001 && occluder.t < source->getDistanceFrom(hit.position)) {
                    power = vec3();
                }
                color_sum += hit.material->shade(hit.position, hit.normal, -ray.dir, source->getLightDirAt(hit.position), power);

            }
            return color_sum;
        }
        return background.getColor(ray.dir);
    }
};

Scene scene;




class FrameBuffer {
    unsigned int textureId;
    vec3 image[windowWidth * windowHeight];

public:
    FrameBuffer() {
        for(int i = 0; i < windowWidth * windowHeight; i++) image[i] = vec3(0.0, 0.0, 0.0);

        glGenTextures(1, &textureId);
        glBindTexture(GL_TEXTURE_2D, textureId);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

    void Bind(Shader* s)
    {
        s->UploadSamplerID();
        glBindTexture(GL_TEXTURE_2D, textureId);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);
    }

    bool ComputeImage()
    {
        static unsigned int iPart = 0;

        if(iPart >= 64)
            return false;
        for(int j = iPart; j < windowHeight; j+=64)
        {
            for(int i = 0; i < windowWidth; i++)
            {
                float ndcX = (2.0 * i - windowWidth) / windowWidth;
                float ndcY = (2.0 * j - windowHeight) / windowHeight;
                Camera& camera = scene.getCamera();
                Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcX, ndcY));

                image[j*windowWidth + i] = scene.trace(ray);
            }
        }
        iPart++;
        return true;
    }
};

class Screen {
    FrameBuffer frameBuffer;
    unsigned int vao;

public:
    Screen()
    {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        unsigned int vbo[2];
        glGenBuffers(2, &vbo[0]);

        glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
        static float vertexCoords[] = { -1, -1,		1, -1,		-1, 1,		1, 1 };

        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);

        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        static float vertexTextureCoords[] = { 0, 0,	1, 0,		0, 1,		1, 1 };

        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexTextureCoords), vertexTextureCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    }

    void Draw(Shader* s)
    {
        if(frameBuffer.ComputeImage())
            glutPostRedisplay();

        s->Run();
        frameBuffer.Bind(s);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisable(GL_BLEND);
    }
};

Screen *screen = 0;


void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    screen->Draw(shader);

    glutSwapBuffers();
}

void onInitialization()
{
    glViewport(0, 0, windowWidth, windowHeight);

    shader = new Shader();

    screen = new Screen();
}

void onExit()
{
    delete screen; screen = 0;
    delete shader; shader = 0;
    printf("exit");
}

int main(int argc, char * argv[]) {
    glutInit(&argc, argv);
#if !defined(__APPLE__)
    glutInitContextVersion(majorVersion, minorVersion);
#endif
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(100, 100);
#if defined(__APPLE__)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE);
#else
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
    glutCreateWindow("Ray Casting");

#if !defined(__APPLE__)
    glewExperimental = true;
    glewInit();
#endif

    printf("GL Vendor    : %s\n", glGetString(GL_VENDOR));
    printf("GL Renderer  : %s\n", glGetString(GL_RENDERER));
    printf("GL Version (string)  : %s\n", glGetString(GL_VERSION));
    glGetIntegerv(GL_MAJOR_VERSION, &majorVersion);
    glGetIntegerv(GL_MINOR_VERSION, &minorVersion);
    printf("GL Version (integer) : %d.%d\n", majorVersion, minorVersion);
    printf("GLSL Version : %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

    glViewport(0, 0, windowWidth, windowHeight);

    onInitialization();

    glutDisplayFunc(onDisplay);

    glutMainLoop();

    onExit();

    return 1;
}



//        for (int i = 0; i < 2; i++) {
//            materials.push_back(new Material(vec3(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX)));
//            vec3 pos = vec3((rand()%9-4)/2.0, (rand()%9-4)/2.0, (rand()%9-4)/6.0);
//            objects.push_back((new Quadric(materials[materials.size()-1]))->ellipsoid()->transform(mat4x4::translation(pos)*mat4x4::scaling(vec3(0.5,0.5,0.5))));
//        }
