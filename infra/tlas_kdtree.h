#pragma once
#include "blas_kdtree.h"

namespace Tmpl8
{
    struct TLASKDTreeNode
    {
        float3 aabbMin = float3(0);
        uint leftRight = 0; // 2x16 bits
        float3 aabbMax = float3(0);
        uint BLAS = 0;
        bool isLeaf() { return leftRight == 0; }
    };

    class TLASKDTree
    {
    private:
        float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
        int FindBestMatch(int* list, int N, int A);
    public:
        TLASKDTree() = default;
        TLASKDTree(std::vector<BLASKDTree*> blasList);
        void Build();
        void Intersect(Ray& ray);
    private:
        TLASKDTreeNode* tlasNode;
        uint nodesUsed = 0, blasCount;
    public:
        std::vector<BLASKDTree*> blas;
        std::chrono::microseconds buildTime;
    };
}
