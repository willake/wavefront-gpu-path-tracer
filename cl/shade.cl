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
  float3 normal;
  float2 uv;
  int matIdx;
} HitInfo;

float3 getFloorNormal() { return (float3)(0, 1, 0); }
float2 getFloorUV(const float3 I) {
  float u = I.x;
  float v = I.z;

  const float invto = 1.0f / 51.2f;

  u *= invto;
  v *= invto;

  u = u - floor(u);
  v = v - floor(v);

  // Return the UV coordinates
  return (float2)(u, v);
}

uint getSkyColor(Ray *ray, uint *pixels, uint width, uint height) {
  if (!pixels)
    return 0;

  float phi = atan2(-ray->D.z, ray->D.x);
  uint u =
      (uint)(width * (phi > 0 ? phi : (phi + 2 * M_PI_F)) * M_2_PI_F - 0.5f);
  uint v = (uint)(height * acos(ray->D.y) * M_1_PI_F - 0.5f);
  uint skyIdx = (u + v * width) % (width * height);

  return pixels[skyIdx];
}

HitInfo getHitInfo(const Ray *ray, const float3 I) {
  HitInfo hitInfo;
  if (ray->objIdx == 1) {
    hitInfo.normal = getFloorNormal();
    hitInfo.uv = getFloorUV(I);
  }
}

__kernel void shade(__global uint *accumulator, __global Ray *rayBuffer,
                    __global uint *skydomePixels, uint skydomeWidth,
                    uint skydomeHeight) {
  // get ray id
  const int index = get_global_id(0);

  Ray ray = rayBuffer[index];

  float3 I = ray.O + ray.D * ray.t;
  HitInfo hitInfo = getHitInfo(&ray, I);

  if (ray.objIdx == -1) {
    accumulator[index] =
        getSkyColor(&ray, skydomePixels, skydomeWidth, skydomeHeight);
  }

  if (ray.objIdx == 1) {
    accumulator[index] = 9527;
  }
}