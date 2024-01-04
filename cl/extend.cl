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
  float v0x, v0y, v0z;
  float v1x, v1y, v1z;
  float v2x, v2y, v2z;
  float cx, cy, cz;
} Tri; // total 48 bytes

typedef struct {
  float n0x, n0y, n0z;
  float n1x, n1y, n1z;
  float n2x, n2y, n2z;
  float uv0x, uv0y, uv0z;
  float uv1x, uv1y, uv1z;
  float uv2x, uv2y, uv2z;
  int dummy;
} TriEx; // total 64 bytes

typedef struct {
  int meshIdx, triStartIdx, triCount;
} MeshInstance; // 12 bytes

typedef struct {
  float aabbMinx, aabbMiny, aabbMinz;
  float aabbMaxx, aabbMaxy, aabbMaxz;
  uint leftFirst, triCount;
} BVHNode; // 32 bytes

typedef struct {
  int meshIdx, startNodeIdx;
  uint nodeCount;
} BVH; // 12 bytes

typedef struct {
  uint objIdx, matIdx, bvhIdx;
  float16 T, invT;
  float aabbMinx, aabbMiny, aabbMinz;
  float aabbMaxx, aabbMaxy, aabbMaxz;
} BLAS;

typedef struct {
  float aabbMinx, aabbMiny, aabbMinz;
  float aabbMaxx, aabbMaxy, aabbMaxz;
  uint leftRight, BLAS;
} TLASNode;

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