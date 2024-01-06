﻿// random numbers: seed using WangHash((threadidx+1)*17), then use RandomInt / RandomFloat
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
    int traversed;      // 4 bytes
    int tested;         // 4 bytes
    bool inside;        // 1 bytes
} Ray;                  // total 77 bytes

typedef struct
{
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

__kernel void testRayStructSize(__global Test *test)
{
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
__kernel void generatePrimaryRays(__global Ray *rayBuffer, __global uint *seeds, __global float4 pixels, int width,
                                  int height, float3 camPos, float3 topLeft, float3 topRight, float3 bottomLeft,
                                  int spp)
{
    // get ray id
    const int index = get_global_id(0);

    const int y = index / width; // Integer division
    const int x = index % width; // Modulo operation

    uint seed = InitSeed(x + y * width + spp * 1799);

    seeds[index] = seed;

    const float u = ((float)x + RandomFloat(seed)) * (1.0f / width);
    const float v = ((float)y + RandomFloat(seed)) * (1.0f / height);
    const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);

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
    pixels[index] = (float4)(1);
}

__kernel void connect()
{
}

__kernel void finalize()
{
}