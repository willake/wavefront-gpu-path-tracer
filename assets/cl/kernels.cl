#include "template/common.h"

// Define ray
typedef struct {
  float3 O, D, rD;    // 36 bytes
  float t;            // 4 bytes
  float2 barycentric; // 8 bytes
  int objIdx;         // 4 bytes
  int triIdx;         // 4 bytes
  int traversed;      // 4 bytes
  int tested;         // 4 bytes
  bool inside;        // 1 bytes
} Ray;                // total 65 bytes

__kernel void generatePrimaryRays(__global Ray *rayBuffer, int width,
                                  int height, float3 camPos, float topLeft,
                                  float topRight, float bottomLeft) {
  // get ray id
  const int index = get_global_id(0);

  const int y = index / width; // Integer division
  const int x = index % width; // Modulo operation

  const float u = (float)x * (1.0f / width);
  const float v = (float)y * (1.0f / height);
  const float3 P =
      topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);

  // initializing a ray
  Ray ray;
  const float3 dir = P - camPos;
  const float distance = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
  ray.O = camPos;
  ray.D = dir / distance;
  ray.rD.x = 1.0f / ray.D.x, ray.rD.y = 1.0f / ray.D.y,
  ray.rD.z = 1.0f / ray.D.z;
  ray.t = 1e34f;
  ray.objIdx = -1, ray.triIdx = -1;

  rayBuffer[index] = ray;
}

__kernel void extend() {}

__kernel void shade() {}

__kernel void connect() {}

__kernel void finalize() {}