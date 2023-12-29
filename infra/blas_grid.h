#pragma once

// reference: https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-acceleration-structure/grid.html
// reference: https://cs184.eecs.berkeley.edu/sp19/lecture/9-44/raytracing
//#define GRID_MAILBOXING // not working very well
namespace Tmpl8
{
	struct GridCell
	{
		std::vector<int> triIndices = {};
	};

	class BLASGrid
	{
	private:
		bool IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		bool IntersectTri(Ray& ray, const Tri& tri, const uint triIdx);
		void IntersectGrid(Ray& ray, long uid);
	public:
		BLASGrid() = default;
		BLASGrid(const int idx, const std::string& modelPath, const mat4 transform, const mat4 scaleMat);
		void Build();
		void Intersect(Ray& ray);
		void SetTransform(mat4 transform);
		float3 GetNormal(const uint triIdx, const float2 barycentric) const;
		float2 GetUV(const uint triIdx, const float2 barycentric) const;
		int GetTriangleCount() const;
	private:
		int3 resolution = 0;
		float3 cellSize = 0;
		std::vector<Tri> triangles;
		std::vector<GridCell> gridCells;
		std::vector<long> mailbox;
		long incremental = 0;
	public:
		int objIdx = -1;
		int matIdx = -1;
		mat4 T, invT;
		aabb localBounds;
		aabb worldBounds;
		std::chrono::microseconds buildTime;
	};
}