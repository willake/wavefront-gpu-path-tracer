#pragma once

#include "base_scene.h"
#include "bvh.h"
#include "grid.h"
#include "model.h"
#include "kdtree.h"
#include "rapidxml.hpp"

//#define USE_BVH
//#define USE_Grid
#define USE_KDTree

#include "tlas_file_scene.h"

namespace Tmpl8
{
	class FileScene : BaseScene
	{
	public:
		FileScene(const string& filePath);
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
#ifdef USE_BVH
		BVH acc;
#endif
#ifdef USE_Grid
		Grid acc;
#endif
#ifdef USE_KDTree
		KDTree acc;
#endif
		std::vector<Model*> models;
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
