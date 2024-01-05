#include "precomp.h"
#include "blas.h"

BLAS::BLAS(int idx, BVH *binded, int materialIdx, mat4 transform)
{
    objIdx = idx;
    bvh = binded;
    matIdx = materialIdx;
    SetTransform(transform);
}

void BLAS::SetTransform(mat4 transform)
{
    T = transform;
    invT = transform.FastInvertedTransformNoScale();
    // update bvh bound
    // calculate world-space bounds using the new matrix
    float3 bmin = bvh->bvhNodes[0].aabbMin, bmax = bvh->bvhNodes[0].aabbMax;
    aabb bounds = aabb();
    for (int i = 0; i < 8; i++)
        bounds.Grow(TransformPosition(float3(i & 1 ? bmax.x : bmin.x, i & 2 ? bmax.y : bmin.y, i & 4 ? bmax.z : bmin.z),
                                      transform));
    aabbMin = bounds.bmin3, aabbMax = bounds.bmax3;
}

float3 BLAS::GetNormal(const uint triIdx, const float2 barycentric) const
{
    return normalize(TransformVector(bvh->GetNormal(triIdx, barycentric), T));
}

float2 BLAS::GetUV(const uint triIdx, const float2 barycentric) const
{
    return bvh->GetUV(triIdx, barycentric);
}

void BLAS::Intersect(Ray &ray)
{
    Ray tRay = Ray(ray);
    tRay.O = TransformPosition(ray.O, invT);
    tRay.D = TransformVector(ray.D, invT);
    tRay.rD = float3(1 / tRay.D.x, 1 / tRay.D.y, 1 / tRay.D.z);

    bvh->IntersectBVH(tRay, objIdx, bvh->rootNodeIdx);

    tRay.O = ray.O;
    tRay.D = ray.D;
    tRay.rD = ray.rD;
    ray = tRay;
}