#pragma once

#include "blas_grid.h"

// reference: https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-acceleration-structure/grid.html
// reference: https://cs184.eecs.berkeley.edu/sp19/lecture/9-44/raytracing
//#define GRID_MAILBOXING // not working very well
namespace Tmpl8
{
	class Grid
	{
	private:
		bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		bool IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectGrid(Ray& ray, long uid);
	public:
		Grid() = default;
		void Build();
		void Intersect(Ray& ray);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
		int3 resolution = 0;
		float3 cellSize = 0;
		std::vector<GridCell> gridCells;
		std::vector<long> mailbox;
		long incremental = 0;
	public:
		int objIdx = -1;
		aabb localBounds;
		std::vector<Tri> triangles;
		std::chrono::microseconds buildTime;
	};
}