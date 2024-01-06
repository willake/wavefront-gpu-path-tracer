#pragma once

#include "base_scene.h"
#define SPEEDTRIX

namespace Tmpl8
{
class PrimitiveScene : BaseScene
{
  protected:
    float3 RandomPointOnLight(const float r0, const float r1) const;

  public:
    PrimitiveScene();
    void SetTime(float t);
    float3 GetSkyColor(const Ray &ray) const;
    float3 GetLightPos() const;
    float3 GetLightColor() const;
    float3 RandomPointOnLight(uint &seed) const;
    void FindNearest(Ray &ray);
    bool IsOccluded(const Ray &ray);
    float3 GetAlbedo(int objIdx, float3 I) const;
    HitInfo GetHitInfo(const Ray &ray, const float3 I);
    int GetTriangleCount() const;

  public:
    Quad quad, dummyQuad1, dummyQuad2, dummyQuad3;
    Sphere sphere;
    Sphere sphere2;
    Cube cube;
    Plane plane[6];
    Torus torus;
    Material materials[11];
};
} // namespace Tmpl8
