#pragma once

#include "blas_bvh.h"

namespace Tmpl8
{
	struct TLASBVHNode
	{
		float3 aabbMin = float3(0);
		uint leftRight = 0; // 2x16 bits
		float3 aabbMax = float3(0);
		uint BLAS = 0;
		bool isLeaf() { return leftRight == 0; }
	};

	class TLASBVH
	{
	private:
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		int FindBestMatch(int* list, int N, int A);
	public:
		TLASBVH() = default;
		TLASBVH(std::vector<BVH*> bvhList);
		void Build();
		void Intersect(Ray& ray);
	private:
		TLASBVHNode* tlasNode;
		uint nodesUsed = 0, blasCount;
	public:
		std::vector<BVH*> blas;
		std::chrono::microseconds buildTime;
	};
}
