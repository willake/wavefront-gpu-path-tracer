#pragma once

#include "blas_grid.h"

namespace Tmpl8
{
    struct TLASGridNode
    {
        float3 aabbMin = float3(0);
        uint leftRight = 0; // 2x16 bits
        float3 aabbMax = float3(0);
        uint BLAS = 0;
        bool isLeaf() { return leftRight == 0; }
    };

    class TLASGrid
    {
    private:
        float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
        int FindBestMatch(int* list, int N, int A);
    public:
        TLASGrid() = default;
        TLASGrid(std::vector<BLASGrid*> blasList);
        void Build();
        void Intersect(Ray& ray);
    private:
        TLASGridNode* tlasNode;
        uint nodesUsed = 0, blasCount;
    public:
        std::vector<BLASGrid*> blas;
        std::chrono::microseconds buildTime;
    };
}
