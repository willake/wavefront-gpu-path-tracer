#pragma once

#define BLAS_BVH_BINS 8

namespace Tmpl8
{
struct GPUBLAS
{
    GPUBLAS(){};
    GPUBLAS(int idx, int bindedBVHIdx, int materialIdx, mat4 transform, float3 boundMin, float3 boundMax)
        : objIdx(idx), matIdx(materialIdx), bvhIdx(bindedBVHIdx), aabbMin(boundMin), aabbMax(boundMax)
    {
        T = transform;
        invT = transform.FastInvertedTransformNoScale();
    };
    uint objIdx = -1, matIdx = -1, bvhIdx = -1;      // 12 bytes
    mat4 T = mat4::Identity();                       // 64 bytes
    mat4 invT = mat4::Identity();                    // 64 bytes
    float3 aabbMin = float3(0), aabbMax = float3(0); // 24 bytes
    // 164 bytes in total
};

// bvh instance
class BLAS
{
  public:
    BLAS(){};
    BLAS(int objIdx, BVH *binded, int matIdx, mat4 transform);
    void Intersect(Ray &ray);
    void SetTransform(mat4 transform);
    float3 GetNormal(const uint triIdx, const float2 barycentric) const;
    float2 GetUV(const uint triIdx, const float2 barycentric) const;

  public:
    int objIdx = -1;
    int matIdx = -1;
    BVH *bvh = nullptr;
    mat4 T = mat4::Identity(), invT = mat4::Identity();
    float3 aabbMin = float3(0), aabbMax = float3(0);
};
} // namespace Tmpl8