#include "precomp.h"
#include "primitive_scene.h"

#define PLANE_X(o, i)                                                                                                  \
    {                                                                                                                  \
        float t = -(ray.O.x + o) * ray.rD.x;                                                                           \
        if (t < ray.t && t > 0)                                                                                        \
            ray.t = t, ray.objIdx = i;                                                                                 \
    }
#define PLANE_Y(o, i)                                                                                                  \
    {                                                                                                                  \
        float t = -(ray.O.y + o) * ray.rD.y;                                                                           \
        if (t < ray.t && t > 0)                                                                                        \
            ray.t = t, ray.objIdx = i;                                                                                 \
    }
#define PLANE_Z(o, i)                                                                                                  \
    {                                                                                                                  \
        float t = -(ray.O.z + o) * ray.rD.z;                                                                           \
        if (t < ray.t && t > 0)                                                                                        \
            ray.t = t, ray.objIdx = i;                                                                                 \
    }

float3 PrimitiveScene::RandomPointOnLight(const float r0, const float r1, uint &lightIdx) const
{
    return float3();
}

PrimitiveScene::PrimitiveScene()
{
    errorMaterial.albedo = float3(255, 192, 203) / 255.f;
    // we store all primitives in one continuous buffer
#ifdef FOURLIGHTS
    for (int i = 0; i < 4; i++)
        quad[i] = Quad(0, 0.5f); // 0: four light sources
#else
    quad = Quad(0, 1); // 0: light source
#endif
    sphere = Sphere(1, float3(0), 0.6f);             // 1: bouncing ball
    sphere2 = Sphere(2, float3(0, 2.5f, -3.07f), 8); // 2: rounded corners
    cube = Cube(3, float3(0), float3(1.15f));        // 3: cube
    plane[0] = Plane(4, float3(1, 0, 0), 3);         // 4: left wall
    plane[1] = Plane(5, float3(-1, 0, 0), 2.99f);    // 5: right wall
    plane[2] = Plane(6, float3(0, 1, 0), 1);         // 6: floor
    plane[3] = Plane(7, float3(0, -1, 0), 2);        // 7: ceiling
    plane[4] = Plane(8, float3(0, 0, 1), 3);         // 8: front wall
    plane[5] = Plane(9, float3(0, 0, -1), 3.99f);    // 9: back wall
    torus = Torus(10, 0.8f, 0.25f);                  // 10: torus
    torus.T = mat4::Translate(-0.25f, 0, 2) * mat4::RotateX(PI / 4);
    torus.invT = torus.T.Inverted();
    materials[0].isLight = true;      // 0: light source
    materials[1].reflectivity = 1.0f; // 1: bouncing ball
    materials[2] = Material();        // 2: rounded corners
    materials[3].refractivity = 1.0f; // 3: cube
    materials[3].absorption = float3(0.5f, 0, 0.5f);
    materials[4] = Material(true); // 4: left wall
    materials[5] = Material(true); // 5: right wall
    materials[6] = Material(true); // 6: floor
    materials[6].reflectivity = 0.3f;
    materials[7] = Material();         // 7: ceiling
    materials[8] = Material();         // 8: front wall
    materials[9] = Material();         // 9: back wall
    materials[10].refractivity = 1.0f; // 10: torus
    SetTime(0);
    // Note: once we have triangle support we should get rid of the class
    // hierarchy: virtuals reduce performance somewhat.
}

void PrimitiveScene::SetTime(float t)
{
    // default time for the scene is simply 0. Updating/ the time per frame
    // enables animation. Updating it per ray can be used for motion blur.
    animTime = t;
#ifdef FOURLIGHTS
    // four light sources are stationary
    quad[0].T = mat4::Translate(-1, 1.5f, -1), quad[0].invT = quad[0].T.FastInvertedTransformNoScale();
    quad[1].T = mat4::Translate(1, 1.5f, -1), quad[1].invT = quad[1].T.FastInvertedTransformNoScale();
    quad[2].T = mat4::Translate(1, 1.5f, 1), quad[2].invT = quad[2].T.FastInvertedTransformNoScale();
    quad[3].T = mat4::Translate(-1, 1.5f, 1), quad[3].invT = quad[3].T.FastInvertedTransformNoScale();
#else
    // light source animation: swing
    mat4 M1base = mat4::Translate(float3(0, 2.6f, 2));
    mat4 M1 = M1base * mat4::RotateZ(sinf(animTime * 0.6f) * 0.1f) * mat4::Translate(float3(0, -0.9f, 0));
    quad.T = M1, quad.invT = M1.FastInvertedTransformNoScale();
#endif
    // cube animation: spin
    mat4 M2base = mat4::RotateX(PI / 4) * mat4::RotateZ(PI / 4);
    mat4 M2 = mat4::Translate(float3(1.8f, 0, 2.5f)) * mat4::RotateY(animTime * 0.5f) * M2base;
    cube.M = M2, cube.invM = M2.FastInvertedTransformNoScale();
    // sphere animation: bounce
    float tm = 1 - sqrf(fmodf(animTime, 2.0f) - 1);
    sphere.pos = float3(-1.8f, -0.4f + tm, 1);
}

float3 PrimitiveScene::GetLightPos() const
{
#ifndef FOURLIGHTS
    // light point position is the middle of the swinging quad
    float3 corner1 = TransformPosition(float3(-0.5f, 0, -0.5f), quad.T);
    float3 corner2 = TransformPosition(float3(0.5f, 0, 0.5f), quad.T);
    return (corner1 + corner2) * 0.5f - float3(0, 0.01f, 0);
#else
    // function is not valid when using four lights; we'll return the origin
    return float3(0);
#endif
}

float3 PrimitiveScene::GetSkyColor(const Ray &ray) const
{
    return float3(0);
}

Light Tmpl8::PrimitiveScene::GetLight(int idx)
{
    return Light();
}

float3 PrimitiveScene::GetLightColor() const
{
    return float3(24, 24, 22);
}

float3 PrimitiveScene::RandomPointOnLight(uint &seed, uint &lightIdx) const
{
    return float3();
}

void PrimitiveScene::FindNearest(Ray &ray)
{
    // room walls - ugly shortcut for more speed
    if (ray.D.x < 0)
        PLANE_X(3, 4) else PLANE_X(-2.99f, 5);
    if (ray.D.y < 0)
        PLANE_Y(1, 6) else PLANE_Y(-2, 7);
    if (ray.D.z < 0)
        PLANE_Z(3, 8) else PLANE_Z(-3.99f, 9);
    quad.Intersect(ray);
    sphere.Intersect(ray);
    sphere2.Intersect(ray);
    cube.Intersect(ray);
    torus.Intersect(ray);
}

bool PrimitiveScene::IsOccluded(const Ray &ray)
{
    if (cube.IsOccluded(ray))
        return true;
    if (sphere.IsOccluded(ray))
        return true;
#ifdef FOURLIGHTS
    for (int i = 0; i < 4; i++)
        if (quad[i].IsOccluded(ray))
            return true;
#else
    if (quad.IsOccluded(ray))
        return true;
#endif
    if (torus.IsOccluded(ray))
        return true;
    return false; // skip planes and rounded corners
}

HitInfo PrimitiveScene::GetHitInfo(const Ray &ray, const float3 I)
{
    HitInfo hitInfo = HitInfo(float3(0), float2(0), &errorMaterial);
    switch (ray.objIdx)
    {
    case 0:
        hitInfo.normal = quad.GetNormal(I);
        hitInfo.uv = float2(0);
        break;
    case 1:
        hitInfo.normal = sphere.GetNormal(I);
        hitInfo.uv = float2(0);
        break;
    case 2:
        hitInfo.normal = sphere2.GetNormal(I);
        hitInfo.uv = float2(0);
        break;
    case 3:
        hitInfo.normal = cube.GetNormal(I);
        hitInfo.uv = float2(0);
        break;
    case 10:
        hitInfo.normal = torus.GetNormal(I);
        hitInfo.uv = float2(0);
        break;
    default:
        hitInfo.normal = float3(0);
        hitInfo.normal[(ray.objIdx - 4) / 2] = 1 - 2 * (float)(ray.objIdx & 1);
        break;
    }
    if (dot(hitInfo.normal, ray.D) > 0)
        hitInfo.normal = -hitInfo.normal;
    hitInfo.material = &materials[ray.objIdx];
    return hitInfo;
}

float3 PrimitiveScene::GetAlbedo(int objIdx, float3 I) const
{
    if (objIdx == -1)
        return float3(0); // or perhaps we should just crash
#ifdef FOURLIGHTS
    if (objIdx == 0)
        return quad[0].GetAlbedo(I); // they're all the same
#else
    if (objIdx == 0)
        return quad.GetAlbedo(I);
#endif
    if (objIdx == 1)
        return sphere.GetAlbedo(I);
    if (objIdx == 2)
        return sphere2.GetAlbedo(I);
    if (objIdx == 3)
        return cube.GetAlbedo(I);
    if (objIdx == 10)
        return torus.GetAlbedo(I);
    return plane[objIdx - 4].GetAlbedo(I);
    // once we have triangle support, we should pass objIdx and the bary-
    // centric coordinates of the hit, instead of the intersection location.
}

int PrimitiveScene::GetTriangleCount() const
{
    return 0;
}