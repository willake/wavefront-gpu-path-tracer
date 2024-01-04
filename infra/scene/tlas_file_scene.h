#pragma once

#include "base_scene.h"
#include "blas_bvh.h"
#include "blas.h"
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
		int meshIdx;
		int materialIdx;
		float3 position;
		float3 rotation;
	};

	struct MeshData {
		std::string modelLocation;
		float3 scale;
	};

	// Define a structure to hold scene information
	struct SceneData {
		std::string name;
		float3 lightPos;
		std::string planeTextureLocation;
		std::string skydomeLocation;
		std::vector<ObjectData> objects;
		std::vector<MeshData> meshes;
		std::vector<MaterialData> materials;
	};


	class TLASFileScene : BaseScene
	{
	public:
		TLASFileScene(const string& filePath);
		SceneData LoadSceneFile(const string& filePath);
		void PrepareBuffers();
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
		TLAS tlas;
		string sceneName;
		Texture skydome;
		Plane floor;
		Quad light;
		int objIdUsed = 2;
		int objCount = 0;
		int materialCount = 0;
		int meshCount = 0;
		int totalTriangleCount = 0;
		int totalBVHNodeCount = 0;
		Material errorMaterial;
		Material primitiveMaterials[3];
		std::vector<Mesh> meshes;
		BVH* bvhs;
		BLAS* blases;
		Material* materials;
		// GPU buffers
		Buffer* skydomeBuffer;
		TextureInfo skydomeInfo;
		Buffer* skydomeInfoBuffer;
		Tri* triangles;
		Buffer* triBuffer;
		TriEx* triangleExs;
		Buffer* triExBuffer;
		uint* triangleIndices;
		Buffer* triIdxBuffer;
		MeshInstance* meshInstances;
		Buffer* meshInsBuffer;
		BVHNode* bvhNodes;
		Buffer* bvhNodeBuffer;
		GPUBVH* gpubvhs;
		Buffer* bvhBuffer;
		GPUBLAS* gpublases;
		Buffer* blasBuffer;
		Buffer* tlasNodeBuffer;
	};
}
