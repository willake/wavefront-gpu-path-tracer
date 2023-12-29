#pragma once

#include "base_scene.h"
#include "blas_bvh.h"
#include "tlas_bvh.h"
#include "rapidxml.hpp"

#define TLAS_USE_BVH

namespace Tmpl8
{
	struct MaterialData {
		float reflectivity;
		float refractivity;
		float3 absorption;
		std::string textureLocation;
	};
	struct ObjectData {
		std::string modelLocation;
		int materialIdx;
		float3 position;
		float3 rotation;
		float3 scale;
	};

	// Define a structure to hold scene information
	struct SceneData {
		std::string name;
		float3 lightPos;
		std::string planeTextureLocation;
		std::string skydomeLocation;
		std::vector<ObjectData> objects;
		std::vector<MaterialData> materials;
	};


	class TLASFileScene : BaseScene
	{
	public:
		TLASFileScene(const string& filePath);
		SceneData LoadSceneFile(const string& filePath);
		void SetTime(float t);
		float3 GetSkyColor(const Ray& ray) const;
		float3 GetLightPos() const;
		float3 GetLightColor() const;
		void FindNearest(Ray& ray);
		bool IsOccluded(const Ray& ray);
		float3 GetAlbedo(int objIdx, float3 I) const;
		HitInfo GetHitInfo(const Ray& ray, const float3 I);
		int GetTriangleCount() const;
		std::chrono::microseconds GetBuildTime() const;
		uint GetMaxTreeDepth() const;
	public:
		float animTime = 0;
		TLASBVH tlas;
		string sceneName;
		Texture skydome;
		Plane floor;
		Quad light;
		int objIdUsed = 2;
		int objCount = 0;
		int materialCount = 0;
		Material errorMaterial;
		Material primitiveMaterials[3];
		std::vector<Material*> materials;
	};
}
