#pragma once

#include "base_scene.h"
#include "blas_bvh.h"
#include "blas.h"
#include "tlas_bvh.h"
#include "rapidxml.hpp"
#include "scene_buffer.h"

#define TLAS_USE_BVH

namespace Tmpl8
{

struct LightData
{
    float3 position;
    float3 rotation;
    float3 color;
    float size;
    // TODO: type
};

struct MaterialData
{
    float reflectivity;
    float refractivity;
    float3 absorption;
    std::string textureLocation;
};

struct ObjectData
{
    int meshIdx;
    int materialIdx;
    float3 position;
    float3 rotation;
};

struct MeshData
{
    std::string modelLocation;
    float3 scale;
};

// Define a structure to hold scene information
struct SceneData
{
    std::string name;
    std::string planeTextureLocation;
    std::string skydomeLocation;
    std::vector<LightData> lights;
    std::vector<ObjectData> objects;
    std::vector<MeshData> meshes;
    std::vector<MaterialData> materials;
};

class TLASFileScene : BaseScene
{
  protected:
    float3 RandomPointOnLight(const float r0, const float r1, uint &lightIdxx) const;

  public:
    TLASFileScene(const string &filePath, SceneBuffer *sceneBuffer = nullptr);
    SceneData LoadSceneFile(const string &filePath);
    void SetTime(float t);
    float3 GetSkyColor(const Ray &ray) const;
    Light GetLightByLightIdx(int lightIdx);
    Light GetLightByObjIdx(int objIdx);
    float3 GetLightPos() const;
    float3 GetLightColor() const;
    float3 RandomPointOnLight(uint &seed, uint &lightIdx) const;
    void FindNearest(Ray &ray);
    bool IsOccluded(const Ray &ray);
    float3 GetAlbedo(int objIdx, float3 I) const;
    HitInfo GetHitInfo(const Ray &ray, const float3 I);
    int GetTriangleCount() const;
    std::chrono::microseconds GetBuildTime() const;
    uint GetMaxTreeDepth() const;

  public:
    float animTime = 0;
    TLAS tlas;
    string sceneName;
    Texture skydome;
    Plane floor;
    Quad *lightQuads;
    Light *lights;
    int objIdUsed = 2;
    uint objCount = 0;
    uint materialCount = 0;
    uint meshCount = 0;
    uint lightCount = 0;
    uint totalTriangleCount = 0;
    uint totalBVHNodeCount = 0;
    uint totalPixelCount = 0;
    Material errorMaterial;
    Material primitiveMaterials[3];
    std::vector<Mesh> meshes;
    BVH *bvhs;
    BLAS *blases;
    Material *materials;
};
} // namespace Tmpl8
