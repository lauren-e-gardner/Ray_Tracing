// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// flag for using soft shadows (set to 1 only when using soft shadows)
#define SOFT_SHADOWS 1

// define number of soft shadow samples to take
#define SOFT_SAMPLING 3.0

// define constant parameters
// EPS is for the precision issue
#define INFINITY 1.0e+12
#define EPS 1.0e-3

// define maximum recursion depth for rays
#define MAX_RECURSION 8

// define constants for scene setting
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
  int shapeType;
  vec3 v1;
  vec3 v2;
  float rad;
};

struct Material {
  int materialType;
  vec3 color;
  float shininess;
  vec3 specular;

  int materialReflectType;
  float reflectivity;
  float refractionRatio;
  int special;
};

struct Object {
  Shape shape;
  Material material;
};

struct Light {
  vec3 position;
  vec3 color;
  float intensity;
  float attenuate;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Intersection {
  vec3 position;
  vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

// varying
varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset(Ray ray, float dist) {
  return ray.origin + (dist * ray.direction);
}

// if a newly found intersection is closer than the best found so far, record
// the new intersection and return true; otherwise leave the best as it was and
// return false.
bool chooseCloserIntersection(float dist, inout float best_dist,
                              inout Intersection intersect,
                              inout Intersection best_intersect) {
  if (best_dist <= dist)
    return false;
  best_dist = dist;
  best_intersect.position = intersect.position;
  best_intersect.normal = intersect.normal;
  return true;
}

// put any general convenience functions you want up here
bool pointInBox(vec2 pmin, vec2 pmax, vec2 pos) {

  bool inbox = false;
  if (lessThanEqual(pos, pmax+EPS)[0] && greaterThanEqual(pos, pmin-EPS)[0]) {
    inbox = true;
  }
  if (lessThanEqual(pos, pmax+EPS)[1] && greaterThanEqual(pos, pmin-EPS)[1] && inbox) {
    inbox = true;
  }
  else {
    inbox = false;
  }
  return inbox;
}
// code derived from: https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
// as linked in assignment website
float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

// forward declaration
float rayIntersectScene(Ray ray, out Material out_mat,
                        out Intersection out_intersect);

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane(Ray ray, vec3 norm, float dist,
                                out Intersection intersect) {
  float a = dot(ray.direction, norm);
  float b = dot(ray.origin, norm) - dist;

  if (a < EPS && a > -EPS)
    return INFINITY;

  float len = -b / a;
  if (len < EPS)
    return INFINITY;

  intersect.position = rayGetOffset(ray, len);
  intersect.normal = norm;
  return len;
}

// Triangle
float findIntersectionWithTriangle(Ray ray, vec3 t1, vec3 t2, vec3 t3,
                                   out Intersection intersect) {

  // intersection with plane variable
  Intersection planeInter;

  // calculating triangle normal, and starting to calculate distance from ray origin to plane
   
  vec3 V = t2 - t1;
  vec3 W = t3 - t1;
  
  vec3 norm = normalize(cross(V, W));
  float d = dot(norm, t1);
  float den = dot(norm, ray.direction);

  // ray is parallel to triangle
  if(den == 0.0) {
    return INFINITY;
  }

  float dist = -(dot(norm, ray.origin) + d)/den; 
  float len = findIntersectionWithPlane(ray, norm, d, planeInter);
  vec3 v3 = t1 - planeInter.position;
  vec3 v4 = t2 - planeInter.position;
  vec3 v5 = t3 - planeInter.position;
  vec3 n1 = cross(v4, v3);
  vec3 n2 = cross(v3, v5);
  vec3 n3 = cross(v5, v4);
  if(dot(ray.direction , n1) < EPS || dot(ray.direction , n2) < EPS || dot(ray.direction , n3) < EPS) {
    return INFINITY;
  }

  else {
    intersect.position =rayGetOffset(ray, len);
    intersect.normal = planeInter.normal;
    return len;
  }
}

// Sphere
float findIntersectionWithSphere(Ray ray, vec3 center, float radius,
                                 out Intersection intersect) {
  
  vec3 L = center - ray.origin;
  float t_ca = dot(L, ray.direction);

  if (t_ca < EPS)
    return INFINITY;

  float d2 = dot(L, L) - (t_ca*t_ca);
  float r2 = (radius*radius);

  if (d2 > r2) {
    return INFINITY;
  }

  float t_hc = sqrt(r2 - d2);
  float t1 = t_ca - t_hc;
  float t2 = t_ca + t_hc;

  if ( t1 > EPS) {
    intersect.position = rayGetOffset(ray, t1);
    intersect.normal = (intersect.position - center) / distance(intersect.position , center);
    return t1;
  }
  else if ( t2 > EPS)  {
    intersect.position = rayGetOffset(ray, t2);
    intersect.normal = (intersect.position - center) / distance(intersect.position , center);
    return t2;
  }
    return INFINITY;
}

// Box
float findIntersectionWithBox(Ray ray, vec3 pmin, vec3 pmax,
                              out Intersection out_intersect) {
  
  // intersection with plane variable
  float close_Intersect = INFINITY;
  Intersection Extra_intersect;

  // normal for each
  vec3 norm;

  //xy
  norm = vec3(0.0, 0.0, 1.0);
  float len = findIntersectionWithPlane(ray, norm, pmin.z, Extra_intersect);
  bool inBox = pointInBox(pmin.xy, pmax.xy, Extra_intersect.position.xy);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.z, Extra_intersect);
  inBox = pointInBox(pmin.xy, pmax.xy, Extra_intersect.position.xy);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }

  //xz
  norm = vec3(0.0, 1.0, 0.0);
  len = findIntersectionWithPlane(ray, norm, pmin.y, Extra_intersect);
  inBox = pointInBox(pmin.xz, pmax.xz, Extra_intersect.position.xz);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.y, Extra_intersect);
  inBox = pointInBox(pmin.xz, pmax.xz, Extra_intersect.position.xz);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }

  //yz
  norm = vec3(1.0, 0.0, 0.0);
  len = findIntersectionWithPlane(ray, norm, pmin.x, Extra_intersect);
  inBox = pointInBox(pmin.yz, pmax.yz, Extra_intersect.position.yz);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.x, Extra_intersect);
  inBox = pointInBox(pmin.yz, pmax.yz, Extra_intersect.position.yz);
  if (len <= close_Intersect && inBox) {
    chooseCloserIntersection(len, close_Intersect, Extra_intersect, out_intersect);
  }
  return close_Intersect;
}

// Cylinder
float getIntersectOpenCylinder(Ray ray, vec3 center, vec3 axis, float len,
                               float rad, out Intersection intersect) {
  return INFINITY;
}

float getIntersectDisc(Ray ray, vec3 center, vec3 norm, float rad,
                       out Intersection intersect) {
  float dist = dot(center, norm);
  Intersection planeIntersect;
  float len = findIntersectionWithPlane(ray, norm, dist, planeIntersect);
  vec3 distFromCenter = planeIntersect.position - center;
  if (dot(norm, distFromCenter) < EPS && pow(length(distFromCenter), 2.0) < (rad * rad)) {
    intersect.position = planeIntersect.position;
    intersect.normal = planeIntersect.normal;
    return len;
  }
  return INFINITY;
}

float findIntersectionWithCylinder(Ray ray, vec3 center, vec3 apex,
                                   float radius,
                                   out Intersection out_intersect) {
  vec3 axis = apex - center;
  float len = length(axis);
  axis = normalize(axis);

  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cylinder
  dist = getIntersectOpenCylinder(ray, center, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- two caps
  dist = getIntersectDisc(ray, center, -axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  dist = getIntersectDisc(ray, apex, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  return best_dist;
}

// Cone
float getIntersectOpenCone(Ray ray, vec3 apex, vec3 axis, float len,
                           float radius, out Intersection intersect) {
  float closest = INFINITY;

  vec3 v_d = ray.origin - apex;
  float phi = dot(v_d, axis);
  float theta = dot(ray.direction, axis);
  float alpha = atan(radius/len);
  float cos_2 = pow(cos(alpha), 2.0);
  float sin_2 = pow(sin(alpha), 2.0);

  float a = pow(length(ray.direction - theta * axis), 2.0) * cos_2 - theta * theta * sin_2;
  float b = 2.0 * (dot(ray.direction - theta * axis, v_d - phi * axis) * cos_2 - theta * phi * sin_2);
  float c = pow(length(v_d - phi * axis), 2.0) * cos_2 - phi * phi * sin_2;

  float quadForm = sqrt(b* b - 4.0 * a * c);

  //quadForm = sqrt(quadForm);
  //positive
  float t = (-b + quadForm)/(2.0 * a);

  vec3 q = rayGetOffset(ray, t);
  vec3 center = apex + axis * len;
  if (dot(axis, q - apex) > EPS && dot(axis, q - center) < -EPS) {
    closest = t;
    intersect.position = q;
  }

// negative
  t = (-b - quadForm)/(2.0*a);
  if (t < EPS) {
    return INFINITY;
  }
  q = rayGetOffset(ray, t);
  if (t < closest && dot(axis, q - apex) > EPS && dot(axis, q - center) < -EPS) {
    closest = t;
    intersect.position = q;
  }


  if ( closest != INFINITY) {
    vec3 apexToPos = intersect.position - apex;
    float x = dot(apexToPos, axis);
    float y = length(apexToPos - x * axis);
    vec3 yhat = normalize(apexToPos - x * axis);
    intersect.normal = normalize(-y * axis + x * yhat);
  }
  return closest;
}

float findIntersectionWithCone(Ray ray, vec3 center, vec3 apex, float radius,
                               out Intersection out_intersect) {
  vec3 axis = center - apex;
  float len = length(axis);
  axis = normalize(axis);

  // -- infinite cone
  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cone
  dist = getIntersectOpenCone(ray, apex, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- caps
  dist = getIntersectDisc(ray, center, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  return best_dist;
}

vec3 calculateSpecialDiffuseColor(Material mat, vec3 posIntersection,
                                  vec3 normalVector) {
  if (mat.special == CHECKERBOARD) {
  float x = floor(posIntersection.x + EPS);
  float y = floor(posIntersection.y + EPS);
  float z = floor(posIntersection.z + EPS);
  float total = ((x + y + z) + EPS);
  if (mod(total, 2.0) > 1.0) {
    mat.color = mat.color * mat.color;
  }
  }
 
   else if (mat.special == MYSPECIAL) {

    float together = mod((-posIntersection.x*posIntersection.y)/2.0, 2.0);
        vec3 black = vec3(0.0,0.0,0.0);
        vec3 green = vec3(0.0,0.5,0.0);
        if ((together > 1.0 + EPS && together < 2.0 + EPS)) mat.color = mat.color;
        else mat.color =  black;
  }

  return mat.color;
  // ----------- STUDENT CODE END ------------
}

vec3 calculateDiffuseColor(Material mat, vec3 posIntersection,
                           vec3 normalVector) {
  // Special colors
  if (mat.special != NONE) {
    return calculateSpecialDiffuseColor(mat, posIntersection, normalVector);
  }
  return vec3(mat.color);
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
bool pointInShadow(vec3 pos, vec3 lightVec) {
  Material hitMaterial;
  Intersection intersect;
  Ray ray;

  //cast ray from pos to light
  ray.origin = pos;
  ray.direction = normalize(lightVec);
  float dist = rayIntersectScene(ray, hitMaterial, intersect);
  float distI = distance(intersect.position, lightVec);
  //float distL = distance( pos - lightVec, ray.origin);

  if(dist + EPS < distI) {
    return true;
  }
  else {
    return false;
  }
}

// use random sampling to compute a ratio that represents the
// fractional contribution of the light to the position pos.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
float softShadowRatio(vec3 pos, vec3 lightVec) {

  Material hitMaterial;
  Intersection intersect;
  Ray ray;
  ray.origin = pos;

  float sampleSize = SOFT_SAMPLING * SOFT_SAMPLING;
  float pi = 3.14159265;

  float inLight = 0.0;

  for (float i = 0.0; i < SOFT_SAMPLING; i++) {
    for (float j = 0.0; j < SOFT_SAMPLING; j++) {
      vec2 noise = vec2((i + rand(pos.xy + vec2(i,j)))/SOFT_SAMPLING, (j + rand(pos.yz + vec2(i,j)))/SOFT_SAMPLING);
      float theta = noise.x * 2.0 * pi;
      float u = noise.y * 2.0 - 1.0;
      vec3 intersectionPoint = vec3(sqrt(1.0 - u * u) * cos(theta), sqrt(1.0 - u * u) * sin(theta), u);
      vec3 offset = lightVec + intersectionPoint;
      ray.direction = normalize(offset);
      float dist = rayIntersectScene(ray, hitMaterial, intersect);

      if (length(offset) - dist < EPS) {
        inLight++;
      }

    }
  }

  return inLight/sampleSize;
}

vec3 getLightContribution(Light light, Material mat, vec3 posIntersection,
                          vec3 normalVector, vec3 eyeVector, bool phongOnly,
                          vec3 diffuseColor) {



  vec3 lightVector = light.position - posIntersection;


  float ratio = 1.0; // default to 1.0 for hard shadows
  if (SOFT_SHADOWS == 1) {
    // if using soft shadows, call softShadowRatio to determine
    // fractional light contribution
    ratio = softShadowRatio(posIntersection, lightVector);
  }
  else {
    // check if point is in shadow with light vector
    if (pointInShadow(posIntersection, lightVector)) {
      return vec3(0.0, 0.0, 0.0);
    }
  }

  // Slight optimization for soft shadows
  if (ratio < EPS) {
    return vec3(0.0, 0.0, 0.0);
  }


  // normalize the light vector for the computations below
  float distToLight = length(lightVector);
  lightVector /= distToLight;

  if (mat.materialType == PHONGMATERIAL ||
      mat.materialType == LAMBERTMATERIAL) {
    vec3 contribution = vec3(0.0, 0.0, 0.0);

    // get light attenuation
    float attenuation = light.attenuate * distToLight;
    float diffuseIntensity =
        max(0.0, dot(normalVector, lightVector)) * light.intensity;

    // glass and mirror objects have specular highlights but no diffuse lighting
    if (!phongOnly) {
      contribution +=
          diffuseColor * diffuseIntensity * light.color / attenuation;
    }

    if (mat.materialType == PHONGMATERIAL) {
      // Start with just black by default (i.e. no Phong term contribution)
      vec3 phongTerm = vec3(0.0, 0.0, 0.0);
      
      vec3 ks = mat.specular;
      float n = mat.shininess;
      vec3 r = reflect(lightVector, normalVector);
      phongTerm = ks * max(0.0, (pow(dot(eyeVector, r), n))) * light.color/light.attenuate;

      contribution += phongTerm;
    }

    return ratio * contribution;
  } else {
    return ratio * diffuseColor;
  }
}

vec3 calculateColor(Material mat, vec3 posIntersection, vec3 normalVector,
                    vec3 eyeVector, bool phongOnly) {
  // The diffuse color of the material at the point of intersection
  // Needed to compute the color when accounting for the lights in the scene
  vec3 diffuseColor = calculateDiffuseColor(mat, posIntersection, normalVector);

  // color defaults to black when there are no lights
  vec3 outputColor = vec3(0.0, 0.0, 0.0);

  // Loop over the MAX_LIGHTS different lights, taking care not to exceed
  // numLights (GLSL restriction), and accumulate each light's contribution
  // to the point of intersection in the scene.
  for (int i = 0; i <= MAX_LIGHTS; i++) {
     if (i >= numLights) {
      break;
    }
    else {
        outputColor = outputColor +  getLightContribution(lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor);
    }
  }
  // Return diffuseColor by default, so you can see something for now.
  return outputColor;
}

// find reflection or refraction direction (depending on material type)
vec3 calcReflectionVector(Material material, vec3 direction, vec3 normalVector,
                          bool isInsideObj) {
  if (material.materialReflectType == MIRRORREFLECT) {
    return reflect(direction, normalVector);
  }

  float eta =
      (isInsideObj) ? 1.0 / material.refractionRatio : material.refractionRatio;
  
float snell = 1.0 - (eta * eta) * (1.0 - dot(normalVector, direction) * dot(normalVector, direction));
vec3 R = vec3(0.0,0.0,0.0);
if(snell >= 0.0) {
  R = eta * direction - (eta*dot(normalVector, direction) + sqrt(snell)) * normalVector;
}
return R;
}

vec3 traceRay(Ray ray) {
  // Accumulate the final color from tracing this ray into resColor.
  vec3 resColor = vec3(0.0, 0.0, 0.0);

  // Accumulate a weight from tracing this ray through different materials
  // based on their BRDFs. Initially all 1.0s (i.e. scales the initial ray's
  // RGB color by 1.0 across all color channels). This captures the BRDFs
  // of the materials intersected by the ray's journey through the scene.
  vec3 resWeight = vec3(1.0, 1.0, 1.0);

  // Flag for whether the ray is currently inside of an object.
  bool isInsideObj = false;

  // Iteratively trace the ray through the scene up to MAX_RECURSION bounces.
  for (int depth = 0; depth < MAX_RECURSION; depth++) {
   
    Material hitMaterial;
    Intersection intersect;
    float dist = rayIntersectScene(ray, hitMaterial, intersect);
    if ((dist >= INFINITY) || (abs(dist) <= EPS) )
    {
      break;
    }

    // Compute the vector from the ray towards the intersection.
    vec3 posIntersection = intersect.position;
    vec3 normalVector    = intersect.normal;

    vec3 eyeVector = normalize(ray.origin - posIntersection);

    // Determine whether we are inside an object using the dot product
    // with the intersection's normal vector
    if (dot(eyeVector, normalVector) < 0.0) {
        normalVector = -normalVector;
        isInsideObj = true;
    } else {
        isInsideObj = false;
    }

    // Material is reflective if it is either mirror or glass in this assignment
    bool reflective = (hitMaterial.materialReflectType == MIRRORREFLECT ||
                       hitMaterial.materialReflectType == GLASSREFLECT);

    // Compute the color at the intersection point based on its material
    // and the lighting in the scene
    vec3 outputColor = calculateColor(hitMaterial, posIntersection,
      normalVector, eyeVector, reflective);

    // A material has a reflection type (as seen above) and a reflectivity
    // attribute. A reflectivity "equal to zero" indicates that the material
    // is neither reflective nor refractive.

    // If a material is neither reflective nor refractive...
    // (1) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (2) Then break the for loop (i.e. do not trace the ray any further).
    if (!reflective) {
      resColor += outputColor * resWeight;
      break;
    }

    // If the material is reflective or refractive...
    // (1) Use calcReflectionVector to compute the direction of the next
    //     bounce of this ray.
    // (2) Update the ray object with the next starting position and
    //     direction to prepare for the next bounce. You should modify the
    //     ray's origin and direction attributes. Be sure to normalize the
    //     direction vector.
    // (3) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (4) Update the current weight using the material's reflectivity
    //     so that it is the appropriate weight for the next ray's color.
    else if (reflective) {
      vec3 direction = calcReflectionVector(hitMaterial, ray.direction, normalVector, isInsideObj);
      ray.origin = intersect.position;
      ray.direction = normalize(direction);
      resColor += outputColor * resWeight;
      resWeight *= hitMaterial.reflectivity;
    }
  }
  return resColor;
}

void main() {
  float cameraFOV = 0.8;
  vec3 direction = vec3(v_position.x * cameraFOV * width / height,
                        v_position.y * cameraFOV, 1.0);

  Ray ray;
  ray.origin = vec3(uMVMatrix * vec4(camera, 1.0));
  ray.direction = normalize(vec3(uMVMatrix * vec4(direction, 0.0)));

  // trace the ray for this pixel
  vec3 res = traceRay(ray);

  // paint the resulting color into this pixel
  gl_FragColor = vec4(res.x, res.y, res.z, 1.0);
}
