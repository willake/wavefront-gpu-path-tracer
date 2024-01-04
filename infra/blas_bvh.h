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

	struct GPUBVH // for GPU
	{
		GPUBVH() {};
		GPUBVH(int meshId, int startIdx, uint count) : meshIdx(meshId), startNodeIdx(startIdx), nodeCount(count) {};
		int meshIdx = -1; // 4 bytes
		int startNodeIdx = -1; // 4 bytes
		uint nodeCount = 0; // 4 bytes
		// 12 bytes in total
	};

	class BVH
	{
	private:
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint nodeIdx, uint depth);
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		void IntersectTri(Ray& ray, const Tri& tri, const int objIdx, const uint triIdx);
		float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
		float CalculateNodeCost(BVHNode& node);
	public:
		BVH() = default;
		BVH(MeshInstance& meshIns, Tri* tri, TriEx* triExs);
		void Build();
		void IntersectBVH(Ray& ray, const int objIdx, const uint nodeIdx);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
	public:
		BVHNode* bvhNodes;
		int triangleCount;
		Tri* triangles;
		TriEx* triangleExs;
		uint* triangleIndices;
		uint rootNodeIdx = 0, nodesUsed = 1;
		std::chrono::microseconds buildTime;
		uint maxDepth = 0;
	};
}