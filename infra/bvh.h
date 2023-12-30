#pragma once

#include "blas_bvh.h"

#define BVH_SAH
#define BVH_FASTER_RAY
#define BVH_BINS 8

// reference: https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/

namespace Tmpl8
{
	class BVH
	{
	private:
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint nodeIdx, uint depth);
#ifdef BVH_FASTER_RAY
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
#else
		bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
#endif
		void IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectBVH(Ray& ray, const uint nodeIdx);
		float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
		float CalculateNodeCost(BVHNode& node);
	public:
		BVH() = default;
		void Build();
		void Refit();
		void Intersect(Ray& ray);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
	public:
		int objIdx = -1;
		std::vector<BVHNode> bvhNodes;
		std::vector<Tri> triangles;
		std::vector<TriEx> triangleExs;
		std::vector<uint> triangleIndices;
		uint rootNodeIdx = 0, nodesUsed = 1;
		std::chrono::microseconds buildTime;
		uint maxDepth = 0;
	};
}