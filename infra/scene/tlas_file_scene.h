#pragma once

#include "base_scene.h"
#include "blas_bvh.h"
#include "blas_grid.h"
#include "blas_kdtree.h"
#include "tlas_bvh.h"
#include "tlas_grid.h"
#include "tlas_kdtree.h"
#include "rapidxml.hpp"

#define TLAS_USE_BVH
//#define TLAS_USE_Grid
//#define TLAS_USE_KDTree

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
#ifdef TLAS_USE_BVH
		TLASBVH tlas;
#endif
#ifdef TLAS_USE_Grid
		TLASGrid tlas;
#endif
#ifdef TLAS_USE_KDTree
		TLASKDTree tlas;
#endif
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
