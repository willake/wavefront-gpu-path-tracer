#pragma once

#include "texture.h"
#include "material.h"

namespace Tmpl8
{
// -----------------------------------------------------------
// Scene class
// We intersect this. The query is internally forwarded to the
// list of primitives, so that the nearest hit can be returned.
// For this hit (distance, obj id), we can query the normal and
// albedo.
// -----------------------------------------------------------
class BaseScene
{
  protected:
    virtual float3 RandomPointOnLight(const float r0, const float r1, uint &lightIdx) const = 0;

  public:
    BaseScene() = default;
    virtual void SetTime(float t) = 0;
    virtual float3 GetSkyColor(const Ray &ray) const = 0;
    virtual Light GetLight(int idx) = 0;
    virtual float3 GetLightPos() const = 0;
    virtual float3 GetLightColor() const = 0;
    virtual float3 RandomPointOnLight(uint &seed, uint &lightIdx) const = 0;
    virtual void FindNearest(Ray &ray) = 0;
    virtual bool IsOccluded(const Ray &ray) = 0;
    virtual float3 GetAlbedo(int objIdx, float3 I) const = 0;
    virtual HitInfo GetHitInfo(const Ray &ray, const float3 I) = 0;
    virtual int GetTriangleCount() const = 0;

  public:
    float animTime = 0;
    Material errorMaterial;
};
} // namespace Tmpl8
