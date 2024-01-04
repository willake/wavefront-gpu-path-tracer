// #include "template/common.h"

// Define ray
// __attribute__((aligned(64)))
typedef struct __attribute__((aligned(128))) {
  float3 O, D, rD;    // 36 bytes
  float t;            // 4 bytes
  float2 barycentric; // 8 bytes
  int objIdx;         // 4 bytes
  int triIdx;         // 4 bytes
  int traversed;
  int tested;
  bool inside; // 1 bytes
} Ray;         // total 64 bytes

typedef struct {
  int sizeO;
  int sizeD;
  int sizerD;
  int sizet;
  int sizebary;
  int sizeObjIdx;
  int sizeTriIdx;
  int sizeInside;
  int sizePadding;
  int sizeTotal;
} Test;

__kernel void testRayStructSize(__global Test *test) {
  const int index = get_global_id(0);
  Test t;
  t.sizeO = sizeof(float3);
  t.sizeD = sizeof(float3);
  t.sizerD = sizeof(float3);
  t.sizet = sizeof(float);
  t.sizebary = sizeof(float2);
  t.sizeObjIdx = sizeof(int);
  t.sizeTriIdx = sizeof(int);
  t.sizeInside = sizeof(bool);
  t.sizePadding = sizeof(char[7]);
  t.sizeTotal = sizeof(Ray);
  test[index] = t;
}
__kernel void generatePrimaryRays(__global Ray *rayBuffer, int width,
                                  int height, float3 camPos, float3 topLeft,
                                  float3 topRight, float3 bottomLeft) {
  // get ray id
  const int index = get_global_id(0);

  const int y = index / width; // Integer division
  const int x = index % width; // Modulo operation

  const float u = (float)x * (1.0f / width);
  const float v = (float)y * (1.0f / height);
  const float3 P =
      topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);

  // // initializing a ray
  Ray ray;
  const float3 dir = P - camPos;
  const float distance = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
  ray.O = camPos;
  ray.D = dir / distance;
  ray.rD = (float3)(1.0f / ray.D.x, 1.0f / ray.D.y, 1.0f / ray.D.z);
  ray.t = 1e34f;
  ray.barycentric = (float2)(0.0f, 0.0f);
  ray.objIdx = -1, ray.triIdx = -1;
  ray.traversed = 0, ray.tested = 0;
  ray.inside = false;

  rayBuffer[index] = ray;
}

__kernel void connect() {}

__kernel void finalize() {}