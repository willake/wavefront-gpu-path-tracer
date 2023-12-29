#include "precomp.h"
#include "kdtree.h"

void KDTree::Build()
{
    auto startTime = std::chrono::high_resolution_clock::now();
    triangleBounds.resize(triangles.size());
    // populate triangle index array
    std::vector<uint> triIndices;
    triIndices.resize(triangles.size());
    for (int i = 0; i < triangles.size(); i++)
    {
        // setup indices
        triIndices[i] = i;
    }
    UpdateBounds();
    // assign all triangles to root node
    rootNode = new KDTreeNode();
    rootNode->aabbMin = localBounds.bmin3;
    rootNode->aabbMax = localBounds.bmax3;
    rootNode->triIndices = triIndices;
    // subdivide recursively
    Subdivide(rootNode, 0);
    auto endTime = std::chrono::high_resolution_clock::now();
    buildTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
}

void KDTree::UpdateBounds()
{
    aabb b;
    for (uint i = 0; i < triangles.size(); i++)
    {
        Tri& tri = triangles[i];
        aabb triBounds;
        triBounds.Grow(tri.vertex0);
        triBounds.Grow(tri.vertex1);
        triBounds.Grow(tri.vertex2);
        b.Grow(triBounds);
        triangleBounds[i] = triBounds;
    }
    localBounds = b;
}


void KDTree::Subdivide(KDTreeNode* node, int depth)
{
    // terminate recursion
    if (depth >= m_maxBuildDepth) return;
    uint triCount = node->triIndices.size();
    if (triCount <= 2) return;
    if (depth > maxDepth) maxDepth = depth;

    // split plane axis and position
    float3 extent = node->aabbMax - node->aabbMin;
    int axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent[axis]) axis = 2;
    float distance = extent[axis] * 0.5f;
    float splitPos = node->aabbMin[axis] + distance;

    std::vector<uint> leftTriIdxs;
    std::vector<uint> rightTriIdxs;

    for (int i = 0; i < triCount; i++)
    {
        uint idx = node->triIndices[i];
        if (triangleBounds[idx].bmax[axis] < splitPos)
        {
            leftTriIdxs.push_back(idx);
        }
        else if (triangleBounds[idx].bmin[axis] > splitPos - 0.001)
        {
            rightTriIdxs.push_back(idx);
        }
        else
        {
            leftTriIdxs.push_back(idx);
            rightTriIdxs.push_back(idx);
        }
    }

    KDTreeNode* leftChild = new KDTreeNode();
    KDTreeNode* rightChild = new KDTreeNode();
    leftChild->triIndices = leftTriIdxs;
    rightChild->triIndices = rightTriIdxs;
    node->left = leftChild;
    node->right = rightChild;
    nodesUsed++; nodesUsed++;

    node->splitAxis = axis;
    node->splitDistance = distance;

    // update the bounds of nodes
    leftChild->aabbMin = node->aabbMin;
    leftChild->aabbMax = node->aabbMax;
    leftChild->aabbMax[axis] = splitPos;

    rightChild->aabbMin = node->aabbMin;
    rightChild->aabbMax = node->aabbMax;
    rightChild->aabbMin[axis] = splitPos;

    node->triIndices.clear();
    node->isLeaf = false;

    // recurse
    Subdivide(leftChild, depth + 1);
    Subdivide(rightChild, depth + 1);
}

bool KDTree::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax, float& tminOut, float& tmaxOut)
{
    float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
    float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
    float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
    tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
    float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
    tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
    tminOut = tmin;
    tmaxOut = tmax;
    return tmax >= tmin && tmin < ray.t && tmax > 0;
}

void KDTree::IntersectTri(Ray& ray, const Tri& tri, const uint triIdx)
{
    const float3 edge1 = tri.vertex1 - tri.vertex0;
    const float3 edge2 = tri.vertex2 - tri.vertex0;
    const float3 h = cross(ray.D, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray.O - tri.vertex0;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray.D, q);
    if (v < 0 || u + v > 1) return;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f)
    {
        if (t < ray.t) ray.t = min(ray.t, t), ray.objIdx = tri.objIdx, ray.triIdx = triIdx, ray.barycentric = float2(u, v);
    }
}

void KDTree::IntersectKDTree(Ray& ray, KDTreeNode* node)
{
    float tmin, tmax;
    if (node == nullptr) return;
    ray.traversed++;
    if (!IntersectAABB(ray, node->aabbMin, node->aabbMax, tmin, tmax)) return;
    if (node->isLeaf)
    {
        uint triCount = node->triIndices.size();
        for (uint i = 0; i < triCount; i++)
        {
            uint triIdx = node->triIndices[i];
            IntersectTri(ray, triangles[triIdx], triIdx);
            ray.tested++;
        }
        return;
    }

    int axis = node->splitAxis;
    float splitPos = node->aabbMin[node->splitAxis] + node->splitDistance;
    float t = (splitPos - ray.O[axis]) / ray.D[axis];

    if (ray.D[axis] > 0)
    {
        // t <= tmin, only the right node is intersected
        if (t < tmin + 0.001)
        {
            IntersectKDTree(ray, node->right);
        }
        // t >= tmax, only the left node is intersected
        else if (t > tmax - 0.001)
        {
            IntersectKDTree(ray, node->left);
        }
        else
        {
            IntersectKDTree(ray, node->left);
            if (ray.t < t) return;
            IntersectKDTree(ray, node->right);
        }
    }
    else
    {
        // t <= tmin, only the left node is intersected
        if (t < tmin + 0.001)
        {
            IntersectKDTree(ray, node->left);
        }
        // t >= tmax, only the left node is intersected
        else if (t > tmax - 0.001)
        {
            IntersectKDTree(ray, node->right);
        }
        else
        {
            IntersectKDTree(ray, node->right);
            if (ray.t < t) return;
            IntersectKDTree(ray, node->left);
        }
    }
}

int KDTree::GetTriangleCount() const
{
    return triangles.size();
}

void KDTree::Intersect(Ray& ray)
{
    IntersectKDTree(ray, rootNode);
}

float3 KDTree::GetNormal(const uint triIdx, const float2 barycentric) const
{
    float3 n0 = triangles[triIdx].normal0;
    float3 n1 = triangles[triIdx].normal1;
    float3 n2 = triangles[triIdx].normal2;
    float3 N = (1 - barycentric.x - barycentric.y) * n0 + barycentric.x * n1 + barycentric.y * n2;
    return normalize(N);
}

float2 KDTree::GetUV(const uint triIdx, const float2 barycentric) const
{
    float2 uv0 = triangles[triIdx].uv0;
    float2 uv1 = triangles[triIdx].uv1;
    float2 uv2 = triangles[triIdx].uv2;
    return (1 - barycentric.x - barycentric.y) * uv0 + barycentric.x * uv1 + barycentric.y * uv2;
}