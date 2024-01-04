#pragma once

#define BLAS_BVH_BINS 8

// reference: https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/

namespace Tmpl8
{
	struct Bin { aabb bounds; int triCount = 0; };

	struct BVHNode
	{
		float3 aabbMin, aabbMax;     // 24 bytes
		uint leftFirst, triCount;   // 8 bytes; total: 32 bytes
		// If it is 0, leftFirst contains the index of the left child node.
		// Otherwise, it contains the index of the first triangle index.
		bool isLeaf() { return triCount > 0; }
	};

	struct BVHInstance
	{
		uint dummy1, dummy2;
		uint idx; // objIdx
		mat4 T;
		mat4 invT; // inverse transform
		uint dummy[6];
	};

	class BLASBVH
	{
	private:
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint nodeIdx, uint depth);
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		void IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectBVH(Ray& ray, const uint nodeIdx);
		float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
		float CalculateNodeCost(BVHNode& node);
	public:
		BLASBVH() = default;
		BLASBVH(const int idx, std::vector<Tri>* tris, std::vector<TriEx>* triExs, const mat4 transform);
		void Build();
		void Refit();
		void Intersect(Ray& ray);
		void SetTransform(mat4 transform);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
	public:
		int objIdx = -1;
		int meshIdx = -1;
		int matIdx = -1;
		std::vector<BVHNode> bvhNodes;
		std::vector<Tri>* triangles;
		std::vector<TriEx>* triangleExs;
		std::vector<uint> triangleIndices;
		uint rootNodeIdx = 0, nodesUsed = 1;
		aabb worldBounds;
		mat4 T, invT;
		std::chrono::microseconds buildTime;
		uint maxDepth = 0;
	};
}