#include "precomp.h"
#include "blas_bvh.h"

BLASBVH::BLASBVH(const int idx, const std::string& modelPath, const mat4 transform, const mat4 scaleMat)
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
			float3(0));
		TriEx triEx(
			vertices[indices[i]].normal,
			vertices[indices[i + 1]].normal,
			vertices[indices[i + 2]].normal,
			vertices[indices[i]].uv,
			vertices[indices[i + 1]].uv,
			vertices[indices[i + 2]].uv,
			objIdx
		);
		tri.centroid = (tri.vertex0 + tri.vertex1 + tri.vertex2) * 0.3333f;
		triangles.push_back(tri);
		triangleExs.push_back(triEx);
	}

	Build();
	SetTransform(transform);
}

void BLASBVH::Build()
{
	auto startTime = std::chrono::high_resolution_clock::now();
	triangleIndices.resize(triangles.size());
	// populate triangle index array
	for (int i = 0; i < triangles.size(); i++)
	{
		// setup indices
		triangleIndices[i] = i;
	}
	// assign all triangles to root node
	bvhNodes.resize(triangles.size() * 2 - 1);
	BVHNode& root = bvhNodes[rootNodeIdx];
	root.leftFirst = 0;
	root.triCount = triangles.size();
	UpdateNodeBounds(rootNodeIdx);
	// subdivide recursively
	Subdivide(rootNodeIdx, 0);
	auto endTime = std::chrono::high_resolution_clock::now();
	buildTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
}

void BLASBVH::Refit()
{
	for (int i = nodesUsed - 1; i >= 0; i--) if (i != 1)
	{
		BVHNode& node = bvhNodes[i];
		if (node.isLeaf())
		{
			// leaf node: adjust bounds to contained triangles
			UpdateNodeBounds(i);
			continue;
		}
		// interior node: adjust bounds to child node bounds
		BVHNode& leftChild = bvhNodes[node.leftFirst];
		BVHNode& rightChild = bvhNodes[node.leftFirst + 1];
		node.aabbMin = fminf(leftChild.aabbMin, rightChild.aabbMin);
		node.aabbMax = fmaxf(leftChild.aabbMax, rightChild.aabbMax);
	}
}

void BLASBVH::UpdateNodeBounds(uint nodeIdx)
{
	BVHNode& node = bvhNodes[nodeIdx];
	node.aabbMin = float3(1e30f);
	node.aabbMax = float3(-1e30f);
	for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		uint leafTriIdx = triangleIndices[first + i];
		Tri& leafTri = triangles[leafTriIdx];
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex0);
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex1);
		node.aabbMin = fminf(node.aabbMin, leafTri.vertex2);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex0);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex1);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.vertex2);
	}
}

void BLASBVH::Subdivide(uint nodeIdx, uint depth)
{
	// terminate recursion
	BVHNode& node = bvhNodes[nodeIdx];
	if (node.triCount <= 2) return;

#ifdef BLAS_BVH_SAH
	// determine split axis using SAH
	int axis;
	float splitPos;
	float splitCost = FindBestSplitPlane(node, axis, splitPos);

	float nosplitCost = CalculateNodeCost(node);
	if (splitCost >= nosplitCost) return;
#else
	// split plane axis and position
	float3 extent = node.aabbMax - node.aabbMin;
	int axis = 0;
	if (extent.y > extent.x) axis = 1;
	if (extent.z > extent[axis]) axis = 2;
	float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
#endif
	// split the group in two halves
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		if (triangles[triangleIndices[i]].centroid[axis] < splitPos)
			i++;
		else
			swap(triangleIndices[i], triangleIndices[j--]);
	}

	// creating child nodes for each half
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == node.triCount) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;

	bvhNodes[leftChildIdx].leftFirst = node.leftFirst;
	bvhNodes[leftChildIdx].triCount = leftCount;
	bvhNodes[rightChildIdx].leftFirst = i;
	bvhNodes[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	UpdateNodeBounds(leftChildIdx);
	UpdateNodeBounds(rightChildIdx);
	if (depth > maxDepth) maxDepth = depth;
	// recurse
	Subdivide(leftChildIdx, depth + 1);
	Subdivide(rightChildIdx, depth + 1);
}

float BLASBVH::CalculateNodeCost(BVHNode& node)
{
	float3 e = node.aabbMax - node.aabbMin; // extent of the node
	float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
	return node.triCount * surfaceArea;
}

float BLASBVH::FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos)
{
	float bestCost = 1e30f;
	for (int a = 0; a < 3; a++)
	{
		float boundsMin = 1e30f, boundsMax = -1e30f;
		for (int i = 0; i < node.triCount; i++)
		{
			Tri& triangle = triangles[triangleIndices[node.leftFirst + i]];
			boundsMin = min(boundsMin, triangle.centroid[a]);
			boundsMax = max(boundsMax, triangle.centroid[a]);
		}
		if (boundsMin == boundsMax) continue;
		// populate the bins
		Bin bin[BLAS_BVH_BINS];
		float scale = BLAS_BVH_BINS / (boundsMax - boundsMin);
		for (uint i = 0; i < node.triCount; i++)
		{
			Tri& triangle = triangles[triangleIndices[node.leftFirst + i]];
			int binIdx = min(BLAS_BVH_BINS - 1,
				(int)((triangle.centroid[a] - boundsMin) * scale));
			bin[binIdx].triCount++;
			bin[binIdx].bounds.Grow(triangle.vertex0);
			bin[binIdx].bounds.Grow(triangle.vertex1);
			bin[binIdx].bounds.Grow(triangle.vertex2);
		}
		// gather data for the 7 planes between the 8 bins
		float leftArea[BLAS_BVH_BINS - 1], rightArea[BLAS_BVH_BINS - 1];
		int leftCount[BLAS_BVH_BINS - 1], rightCount[BLAS_BVH_BINS - 1];
		aabb leftBox, rightBox;
		int leftSum = 0, rightSum = 0;
		for (int i = 0; i < BLAS_BVH_BINS - 1; i++)
		{
			leftSum += bin[i].triCount;
			leftCount[i] = leftSum;
			leftBox.Grow(bin[i].bounds);
			leftArea[i] = leftBox.Area();
			rightSum += bin[BLAS_BVH_BINS - 1 - i].triCount;
			rightCount[BLAS_BVH_BINS - 2 - i] = rightSum;
			rightBox.Grow(bin[BLAS_BVH_BINS - 1 - i].bounds);
			rightArea[BLAS_BVH_BINS - 2 - i] = rightBox.Area();
		}
		// calculate SAH cost for the 7 planes
		scale = (boundsMax - boundsMin) / BLAS_BVH_BINS;
		for (int i = 0; i < BLAS_BVH_BINS - 1; i++)
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

#ifdef BLAS_BVH_FASTER_RAY
float BLASBVH::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}
#else
bool BLASBVH::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	return tmax >= tmin && tmin < ray.t && tmax > 0;
}
#endif
void BLASBVH::IntersectTri(Ray& ray, const Tri& tri, const uint triIdx)
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

void BLASBVH::IntersectBVH(Ray& ray, const uint nodeIdx)
{
#ifdef BLAS_BVH_FASTER_RAY
	BVHNode* node = &bvhNodes[rootNodeIdx], * stack[64];
	uint stackPtr = 0;
	while (1)
	{
		ray.traversed++;
		if (node->isLeaf())
		{
			for (uint i = 0; i < node->triCount; i++)
			{
				uint triIdx = triangleIndices[node->leftFirst + i];
				ray.tested++;
				IntersectTri(ray, triangles[triIdx], triIdx);
			}
			if (stackPtr == 0) break; else node = stack[--stackPtr];

			continue;
		}
		BVHNode* child1 = &bvhNodes[node->leftFirst];
		BVHNode* child2 = &bvhNodes[node->leftFirst + 1];
		float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
		float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax);
		if (dist1 > dist2) { swap(dist1, dist2); swap(child1, child2); }
		if (dist1 == 1e30f)
		{
			if (stackPtr == 0) break; else node = stack[--stackPtr];
		}
		else
		{
			node = child1;
			if (dist2 != 1e30f) stack[stackPtr++] = child2;
		}
	}
#else
	BVHNode& node = bvhNodes[nodeIdx];
	ray.traversed++;
	if (!IntersectAABB(ray, node.aabbMin, node.aabbMax)) return;
	if (node.isLeaf())
	{
		for (uint i = 0; i < node.triCount; i++)
		{
			uint triIdx = triangleIndices[node.leftFirst + i];
			ray.tested++;
			IntersectTri(ray, triangles[triIdx], triIdx);
		}
	}
	else
	{
		IntersectBVH(ray, node.leftFirst);
		IntersectBVH(ray, node.leftFirst + 1);
	}
#endif
}

int BLASBVH::GetTriangleCount() const
{
	return triangles.size();
}

void BLASBVH::SetTransform(mat4 transform)
{
	T = transform;
	invT = transform.FastInvertedTransformNoScale();
	// update bvh bound
	// calculate world-space bounds using the new matrix
	float3 bmin = bvhNodes[0].aabbMin, bmax = bvhNodes[0].aabbMax;
	worldBounds = aabb();
	for (int i = 0; i < 8; i++)
		worldBounds.Grow(TransformPosition(float3(i & 1 ? bmax.x : bmin.x,
			i & 2 ? bmax.y : bmin.y, i & 4 ? bmax.z : bmin.z), transform));
}

void BLASBVH::Intersect(Ray& ray)
{
	Ray tRay = Ray(ray);
	tRay.O = TransformPosition_SSE(ray.O4, invT);
	tRay.D = TransformVector_SSE(ray.D4, invT);
	tRay.rD = float3(1 / tRay.D.x, 1 / tRay.D.y, 1 / tRay.D.z);

	IntersectBVH(tRay, rootNodeIdx);

	tRay.O = ray.O;
	tRay.D = ray.D;
	tRay.rD = ray.rD;
	ray = tRay;
}

float3 BLASBVH::GetNormal(const uint triIdx, const float2 barycentric) const
{
	float3 n0 = triangleExs[triIdx].normal0;
	float3 n1 = triangleExs[triIdx].normal1;
	float3 n2 = triangleExs[triIdx].normal2;
	float3 N = (1 - barycentric.x - barycentric.y) * n0 + barycentric.x * n1 + barycentric.y * n2;
	return normalize(TransformVector(N, T));
}

float2 BLASBVH::GetUV(const uint triIdx, const float2 barycentric) const
{
	float2 uv0 = triangleExs[triIdx].uv0;
	float2 uv1 = triangleExs[triIdx].uv1;
	float2 uv2 = triangleExs[triIdx].uv2;
	return (1 - barycentric.x - barycentric.y) * uv0 + barycentric.x * uv1 + barycentric.y * uv2;
}