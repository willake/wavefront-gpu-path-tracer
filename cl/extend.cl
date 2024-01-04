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

void intersectFloor(Ray *ray) {
  float3 N = (float3)(0, 1, 0); // dist = 1
  float t = -(dot(ray->O, N) + 1) / (dot(ray->D, N));
  if (t < ray->t && t > 0)
    ray->t = t, ray->objIdx = 1;
}

__kernel void extend(__global Ray *rayBuffer) {
  const int index = get_global_id(0);

  Ray ray = rayBuffer[index];

  intersectFloor(&ray);

  rayBuffer[index] = ray;
}