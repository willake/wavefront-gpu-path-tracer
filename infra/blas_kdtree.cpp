#include "precomp.h"
#include "blas_kdtree.h"

BLASKDTree::BLASKDTree(const int idx, const std::string& modelPath, const mat4 transform, const mat4 scaleMat)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, modelPath.c_str()))
    {
        throw std::runtime_error(warn + err);
    }

    std::unordered_map<Vertex, uint32_t> uniqueVertices{};

    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;

    for (const auto& shape : shapes)
    {
        for (const auto& index : shape.mesh.indices)
        {
            Vertex vertex{};

            if (index.vertex_index >= 0)
            {
                vertex.position = {
                    attrib.vertices[3 * index.vertex_index + 0],
                    attrib.vertices[3 * index.vertex_index + 1],
                    attrib.vertices[3 * index.vertex_index + 2] };
            }

            if (index.normal_index >= 0)
            {
                vertex.normal = {
                    attrib.normals[3 * index.normal_index + 0],
                    attrib.normals[3 * index.normal_index + 1],
                    attrib.normals[3 * index.normal_index + 2] };
            }

            if (index.texcoord_index >= 0)
            {
                vertex.uv = {
                    attrib.texcoords[2 * index.texcoord_index + 0],
                    attrib.texcoords[2 * index.texcoord_index + 1] };
            }

            if (uniqueVertices.count(vertex) == 0)
            {
                uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
                vertices.push_back(vertex);
            }
            indices.push_back(uniqueVertices[vertex]);
        }
    }

    objIdx = idx;

    for (int i = 0; i < indices.size(); i += 3)
    {
        Tri tri(
            TransformPosition(vertices[indices[i]].position, scaleMat),
            TransformPosition(vertices[indices[i + 1]].position, scaleMat),
            TransformPosition(vertices[indices[i + 2]].position, scaleMat),
            vertices[indices[i]].normal,
            vertices[indices[i + 1]].normal,
            vertices[indices[i + 2]].normal,
            vertices[indices[i]].uv,
            vertices[indices[i + 1]].uv,
            vertices[indices[i + 2]].uv,
            float3(0), objIdx);
        tri.centroid = (tri.vertex0 + tri.vertex1 + tri.vertex2) * 0.3333f;
        triangles.push_back(tri);
    }

    Build();
    SetTransform(transform);
}

void BLASKDTree::Build()
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

void BLASKDTree::UpdateBounds()
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

float BLASKDTree::EvaluateSAH(KDTreeNode* node, int axis, float pos)
{
    // determine triangle counts and bounds for this split candidate
    aabb leftBox, rightBox;
    int leftCount = 0, rightCount = 0;
    int triCount = node->triIndices.size();
    for (uint i = 0; i < triCount; i++)
    {
        Tri& triangle = triangles[node->triIndices[i]];
        aabb bounds = triangleBounds[node->triIndices[i]];
        if (bounds.bmax[axis] < pos + 0.001)
        {
            leftCount++;
            leftBox.Grow(triangle.vertex0);
            leftBox.Grow(triangle.vertex1);
            leftBox.Grow(triangle.vertex2);
        }
        else if (bounds.bmin[axis] > pos - 0.001)
        {
            rightCount++;
            rightBox.Grow(triangle.vertex0);
            rightBox.Grow(triangle.vertex1);
            rightBox.Grow(triangle.vertex2);
        }
        else
        {
            leftCount++;
            leftBox.Grow(triangle.vertex0);
            leftBox.Grow(triangle.vertex1);
            leftBox.Grow(triangle.vertex2);
            rightCount++;
            rightBox.Grow(triangle.vertex0);
            rightBox.Grow(triangle.vertex1);
            rightBox.Grow(triangle.vertex2);
        }
    }
    float cost = leftCount * leftBox.Area() + rightCount * rightBox.Area();

    return cost > 0 ? cost : 1e30f;
}

float BLASKDTree::CalculateNodeCost(KDTreeNode* node)
{
    float3 e = node->aabbMax - node->aabbMin; // extent of the node
    float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
    return node->triIndices.size() * surfaceArea;
}

float BLASKDTree::FindBestSplitPlane(KDTreeNode* node, int& axis, float& splitPos)
{
    float bestCost = 1e30f;
    int triCount = node->triIndices.size();
    for (int a = 0; a < 3; a++)
    {
        float boundsMin = 1e30f, boundsMax = -1e30f;
        for (int i = 0; i < triCount; i++)
        {
            Tri& triangle = triangles[node->triIndices[i]];
            boundsMin = min(boundsMin, triangle.centroid[a]);
            boundsMax = max(boundsMax, triangle.centroid[a]);
        }
        if (boundsMin == boundsMax) continue;
        // populate the bins
        Bin bin[KD_BINS];
        float scale = KD_BINS / (boundsMax - boundsMin);
        for (uint i = 0; i < triCount; i++)
        {
            Tri& triangle = triangles[node->triIndices[i]];
            int binIdx = min(KD_BINS - 1,
                (int)((triangle.centroid[a] - boundsMin) * scale));
            bin[binIdx].triCount++;
            bin[binIdx].bounds.Grow(triangle.vertex0);
            bin[binIdx].bounds.Grow(triangle.vertex1);
            bin[binIdx].bounds.Grow(triangle.vertex2);
        }
        // gather data for the 7 planes between the 8 bins
        float leftArea[KD_BINS - 1], rightArea[KD_BINS - 1];
        int leftCount[KD_BINS - 1], rightCount[KD_BINS - 1];
        aabb leftBox, rightBox;
        int leftSum = 0, rightSum = 0;
        for (int i = 0; i < KD_BINS - 1; i++)
        {
            leftSum += bin[i].triCount;
            leftCount[i] = leftSum;
            leftBox.Grow(bin[i].bounds);
            leftArea[i] = leftBox.Area();
            rightSum += bin[KD_BINS - 1 - i].triCount;
            rightCount[KD_BINS - 2 - i] = rightSum;
            rightBox.Grow(bin[KD_BINS - 1 - i].bounds);
            rightArea[KD_BINS - 2 - i] = rightBox.Area();
        }
        // calculate SAH cost for the 7 planes
        scale = (boundsMax - boundsMin) / KD_BINS;
        for (int i = 0; i < KD_BINS - 1; i++)
        {
            float planeCost =
                leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
            if (planeCost < bestCost)
                axis = a, splitPos = boundsMin + scale * (i + 1),
                bestCost = planeCost;
        }
    }
    return bestCost;
}

void BLASKDTree::Subdivide(KDTreeNode* node, int depth)
{
    // terminate recursion
    if (depth >= m_maxBuildDepth) return;
    uint triCount = node->triIndices.size();
    if (triCount <= 2) return;

#ifdef KD_SAH
    // determine split axis using SAH
    int axis;
    float splitPos;
    float splitCost = FindBestSplitPlane(node, axis, splitPos);

    float nosplitCost = CalculateNodeCost(node);
    if (splitCost >= nosplitCost) return;
    float distance = splitPos - node->aabbMin[axis];
#else
    // split plane axis and position
    float3 extent = node->aabbMax - node->aabbMin;
    int axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent[axis]) axis = 2;
    float distance = extent[axis] * 0.5f;
    float splitPos = node->aabbMin[axis] + distance;
#endif
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

    //if (leftTriIdxs.size() == triCount || rightTriIdxs.size() == triCount) return;

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

bool BLASKDTree::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax, float& tminOut, float& tmaxOut)
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

void BLASKDTree::IntersectTri(Ray& ray, const Tri& tri, const uint triIdx)
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
        if (t < ray.t) ray.t = min(ray.t, t), ray.objIdx = objIdx, ray.triIdx = triIdx, ray.barycentric = float2(u, v);
    }
}

void BLASKDTree::IntersectKDTree(Ray& ray, KDTreeNode* node)
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

    //IntersectKDTree(ray, node->left);
    //IntersectKDTree(ray, node->right);

    if(ray.D[axis] > 0)
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
            if (ray.objIdx == objIdx && ray.t < t) return;
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
            if (ray.objIdx == objIdx && ray.t < t) return;
            IntersectKDTree(ray, node->left);
        }
    }
}

int BLASKDTree::GetTriangleCount() const
{
    return triangles.size();
}

void BLASKDTree::SetTransform(mat4 transform)
{
    T = transform;
    invT = transform.FastInvertedTransformNoScale();
    // update bvh bound
    // calculate world-space bounds using the new matrix
    float3 bmin = rootNode->aabbMin, bmax = rootNode->aabbMax;
    worldBounds = aabb();
    for (int i = 0; i < 8; i++)
        worldBounds.Grow(TransformPosition(float3(i & 1 ? bmax.x : bmin.x,
            i & 2 ? bmax.y : bmin.y, i & 4 ? bmax.z : bmin.z), transform));
}

void BLASKDTree::Intersect(Ray& ray)
{
    Ray tRay = Ray(ray);
    tRay.O = TransformPosition_SSE(ray.O4, invT);
    tRay.D = TransformVector_SSE(ray.D4, invT);
    tRay.rD = float3(1 / tRay.D.x, 1 / tRay.D.y, 1 / tRay.D.z);

    IntersectKDTree(tRay, rootNode);

    tRay.O = ray.O;
    tRay.D = ray.D;
    tRay.rD = ray.rD;
    ray = tRay;
}

float3 BLASKDTree::GetNormal(const uint triIdx, const float2 barycentric) const
{
    float3 n0 = triangles[triIdx].normal0;
    float3 n1 = triangles[triIdx].normal1;
    float3 n2 = triangles[triIdx].normal2;
    float3 N = (1 - barycentric.x - barycentric.y) * n0 + barycentric.x * n1 + barycentric.y * n2;
    return normalize(TransformVector(N, T));
}

float2 BLASKDTree::GetUV(const uint triIdx, const float2 barycentric) const
{
    float2 uv0 = triangles[triIdx].uv0;
    float2 uv1 = triangles[triIdx].uv1;
    float2 uv2 = triangles[triIdx].uv2;
    return (1 - barycentric.x - barycentric.y) * uv0 + barycentric.x * uv1 + barycentric.y * uv2;
}