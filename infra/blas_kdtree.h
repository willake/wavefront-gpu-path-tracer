#pragma once

//#define KD_SAH
#define KD_FASTER_RAY
#define KD_BINS 8
// reference: course slides
// reference: https://www.youtube.com/watch?v=TrqK-atFfWY&ab_channel=JustinSolomon
// reference: https://github.com/reddeupenn/kdtreePathTracerOptimization
// refernece: On building fast kd-Trees for Ray Tracing, and on doing that in O(N log N) by Ingo Wald and Vlastimil Havran

// This is the KDTree class specially for TLAS, which a KDTree represent a model

namespace Tmpl8
{
	struct KDTreeNode
	{
		float3 aabbMin = 0, aabbMax = 0;
		KDTreeNode* left;
		KDTreeNode* right;
		int splitAxis = 0;
		float splitDistance = 0;
		std::vector<uint> triIndices;
		bool isLeaf = true;
	};

	class BLASKDTree
	{
	private:
		void UpdateBounds();
		float CalculateNodeCost(KDTreeNode* node);
		float FindBestSplitPlane(KDTreeNode* node, int& axis, float& splitPos);
		void Subdivide(KDTreeNode* node, int depth);
		float EvaluateSAH(KDTreeNode* node, int axis, float pos);
		bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax, float& tmin, float& tmax);
		void IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectKDTree(Ray& ray, KDTreeNode* node);
	public:
		BLASKDTree() = default;
		BLASKDTree(const int idx, const std::string& modelPath, const mat4 transform, const mat4 scaleMat);
		void Build();
		void Intersect(Ray& ray);
		void SetTransform(mat4 transform);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
		int m_maxBuildDepth = 20;
	public:
		int objIdx = -1;
		int matIdx = -1;
		KDTreeNode* rootNode;
		std::vector<Tri> triangles;
		std::vector<aabb> triangleBounds;
		uint rootNodeIdx = 0, nodesUsed = 1;
		aabb localBounds;
		aabb worldBounds;
		mat4 T, invT;
		std::chrono::microseconds buildTime;
	};
}