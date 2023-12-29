#pragma once

#include "blas_kdtree.h"
// reference: course slides
// reference: https://www.youtube.com/watch?v=TrqK-atFfWY&ab_channel=JustinSolomon
// reference: https://github.com/reddeupenn/kdtreePathTracerOptimization
// refernece: On building fast kd-Trees for Ray Tracing, and on doing that in O(N log N) by Ingo Wald and Vlastimil Havran

// This is the KDTree class specially for representing a scene

namespace Tmpl8
{
	class KDTree
	{
	private:
		void UpdateBounds();
		void Subdivide(KDTreeNode* node, int depth);
		bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax, float& tmin, float& tmax);
		void IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectKDTree(Ray& ray, KDTreeNode* node);
	public:
		KDTree() = default;
		void Build();
		void Intersect(Ray& ray);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
		int m_maxBuildDepth = 20;
	public:
		KDTreeNode* rootNode;
		std::vector<Tri> triangles;
		std::vector<aabb> triangleBounds;
		uint rootNodeIdx = 0, nodesUsed = 1;
		aabb localBounds;
		std::chrono::microseconds buildTime;
		uint maxDepth = 0;
	};
}