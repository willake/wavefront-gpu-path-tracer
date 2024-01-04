#pragma once

#include "blas_bvh.h"

namespace Tmpl8
{
	struct TLASNode
	{
		float3 aabbMin = float3(0); // 12 bytes
		uint leftRight = 0; // 4 bytes
		float3 aabbMax = float3(0); // 12 bytes
		uint BLAS = 0; // 4 bytes
		// 32 bytes in total
		bool isLeaf() { return leftRight == 0; }
	};

	class TLAS
	{
	private:
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		int FindBestMatch(int* list, int N, int A);
	public:
		TLAS() = default;
		TLAS(BLAS* blasList, uint count);
		void Build();
		void Intersect(Ray& ray);
	private:
	public:
		TLASNode* tlasNode;
		BLAS* blases;
		uint nodesUsed = 0, blasCount;
		std::chrono::microseconds buildTime;
	};
}
