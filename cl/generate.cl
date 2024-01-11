// random numbers: seed using WangHash((threadidx+1)*17), then use RandomInt / RandomFloat
uint WangHash(uint s)
{
    s = (s ^ 61) ^ (s >> 16), s *= 9, s = s ^ (s >> 4), s *= 0x27d4eb2d, s = s ^ (s >> 15);
    return s;
}
uint RandomInt(uint *s)
{
    *s ^= *s << 13, *s ^= *s >> 17, *s ^= *s << 5;
    return *s;
}
float RandomFloat(uint *s)
{
    return RandomInt(s) * 2.3283064365387e-10f; /* = 1 / (2^32-1) */
}

uint InitSeed(uint seedBase)
{
    return WangHash((seedBase + 1) * 17);
}

// Define ray
// __attribute__((aligned(64)))
typedef struct __attribute__((aligned(128)))
{
    float3 O, D, rD;    // 48 bytes
    float2 barycentric; // 8 bytes
    float t;            // 4 bytes
    int objIdx;         // 4 bytes
    int triIdx;         // 4 bytes
    int pixelIdx;       // 4 bytes
    int traversed;      // 4 bytes
    int tested;         // 4 bytes
    bool inside;        // 1 bytes
    bool lastSpecular;  // 1 bytes
} Ray;                  // total 82 bytes

__kernel void generatePrimaryRays(__global Ray *rayBuffer, __global uint *seeds, __global float4 *pixels, int width,
                                  int height, float3 camPos, float3 topLeft, float3 topRight, float3 bottomLeft,
                                  int spp)
{
    // get ray id
    const int index = get_global_id(0);

    const int y = index / width; // Integer division
    const int x = index % width; // Modulo operation

    uint seed = InitSeed(index + spp * 1799);
    float r = RandomFloat(&seed);
    const float u = ((float)x + r) * (1.0f / width);
    const float v = ((float)y + r) * (1.0f / height);
    const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);

    // // // initializing a ray
    Ray ray;
    const float3 dir = P - camPos;
    ray.O = camPos;
    ray.D = normalize(dir);
    ray.rD = (float3)(1.0f / ray.D.x, 1.0f / ray.D.y, 1.0f / ray.D.z);
    ray.t = 1e34f;
    ray.barycentric = (float2)(0.0f, 0.0f);
    ray.objIdx = -1, ray.triIdx = -1;
    ray.traversed = 0, ray.tested = 0;
    ray.inside = false;
    ray.pixelIdx = index;

    rayBuffer[index] = ray;
    seeds[index] = seed;
    pixels[index] = (float4)(1);
}

__kernel void clearAccumulator(__global float4 *accumulator)
{
    const int index = get_global_id(0);
    accumulator[index] = (float4)(0);
}