#include "precomp.h"
#include "tlas_bvh.h"

TLAS::TLAS(BLAS* blasList, uint count)
{
	blases = blasList;
	blasCount = count;
	// allocate TLAS nodes
	tlasNode = (TLASNode*)_aligned_malloc(sizeof(TLASNode) * 2 * blasCount, 64);
	Build();
}

void TLAS::Build()
{
	auto startTime = std::chrono::high_resolution_clock::now();
	// assign a TLASleaf node to each BLAS
	int nodeIdx[256], nodeIndices = blasCount;
	nodesUsed = 1;
	for (uint i = 0; i < blasCount; i++)
	{
		nodeIdx[i] = nodesUsed;
		tlasNode[nodesUsed].aabbMin = blases[i].worldBounds.bmin3;
		tlasNode[nodesUsed].aabbMax = blases[i].worldBounds.bmax3;
		tlasNode[nodesUsed].BLAS = i;
		tlasNode[nodesUsed++].leftRight = 0; // makes it a leaf
	}

	// use agglomerative clustering to build the TLAS
	int A = 0, B = FindBestMatch(nodeIdx, nodeIndices, A);
	while (nodeIndices > 1)
	{
		int C = FindBestMatch(nodeIdx, nodeIndices, B);
		if (A == C)
		{
			int nodeIdxA = nodeIdx[A], nodeIdxB = nodeIdx[B];
			TLASNode& nodeA = tlasNode[nodeIdxA];
			TLASNode& nodeB = tlasNode[nodeIdxB];
			TLASNode& newNode = tlasNode[nodesUsed];
			newNode.leftRight = nodeIdxA + (nodeIdxB << 16);
			newNode.aabbMin = fminf(nodeA.aabbMin, nodeB.aabbMin);
			newNode.aabbMax = fmaxf(nodeA.aabbMax, nodeB.aabbMax);
			nodeIdx[A] = nodesUsed++;
			nodeIdx[B] = nodeIdx[nodeIndices - 1];
			B = FindBestMatch(nodeIdx, --nodeIndices, A);
		}
		else A = B, B = C;
	}
	tlasNode[0] = tlasNode[nodeIdx[A]];
	auto endTime = std::chrono::high_resolution_clock::now();
	buildTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
}

int TLAS::FindBestMatch(int* list, int N, int A)
{
	float smallest = 1e30f;
	int bestB = -1;
	for (int B = 0; B < N; B++) if (B != A)
	{
		float3 bmax = fmaxf(tlasNode[list[A]].aabbMax, tlasNode[list[B]].aabbMax);
		float3 bmin = fminf(tlasNode[list[A]].aabbMin, tlasNode[list[B]].aabbMin);
		float3 e = bmax - bmin;
		float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
		if (surfaceArea < smallest) smallest = surfaceArea, bestB = B;
	}
	return bestB;
}

float TLAS::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}

void TLAS::Intersect(Ray& ray)
{
	TLASNode* node = &tlasNode[0], * stack[64];
	uint stackPtr = 0;
	while (1)
	{
		ray.traversed++;
		if (node->isLeaf())
		{
			blases[node->BLAS].Intersect(ray);
			if (stackPtr == 0) break; else node = stack[--stackPtr];
			continue;
		}
		TLASNode* child1 = &tlasNode[node->leftRight & 0xffff];
		TLASNode* child2 = &tlasNode[node->leftRight >> 16];
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
}