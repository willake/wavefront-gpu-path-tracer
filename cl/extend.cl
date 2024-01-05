// Define ray
// __attribute__((aligned(64)))
typedef struct __attribute__((aligned(128)))
{
    float3 O, D, rD;    // 36 bytes
    float t;            // 4 bytes
    float2 barycentric; // 8 bytes
    int objIdx;         // 4 bytes
    int triIdx;         // 4 bytes
    int traversed;
    int tested;
    bool inside; // 1 bytes
} Ray;           // total 64 bytes

typedef struct
{
    float v0x, v0y, v0z;
    float v1x, v1y, v1z;
    float v2x, v2y, v2z;
    float cx, cy, cz;
} Tri; // total 48 bytes

typedef struct
{
    float n0x, n0y, n0z;
    float n1x, n1y, n1z;
    float n2x, n2y, n2z;
    float uv0x, uv0y, uv0z;
    float uv1x, uv1y, uv1z;
    float uv2x, uv2y, uv2z;
    int dummy;
} TriEx; // total 64 bytes

typedef struct
{
    int meshIdx, triStartIdx, triCount;
} MeshInstance; // 12 bytes

typedef struct
{
    float aabbMinx, aabbMiny, aabbMinz;
    float aabbMaxx, aabbMaxy, aabbMaxz;
    uint leftFirst, triCount;
} BVHNode; // 32 bytes

typedef struct
{
    int meshIdx, startNodeIdx;
    uint nodeCount;
} BVH; // 12 bytes

typedef struct
{
    uint objIdx, matIdx, bvhIdx;
    float16 T, invT;
    float aabbMinx, aabbMiny, aabbMinz;
    float aabbMaxx, aabbMaxy, aabbMaxz;
} BLAS;

typedef struct
{
    float aabbMinx, aabbMiny, aabbMinz;
    float aabbMaxx, aabbMaxy, aabbMaxz;
    uint leftRight, BLAS;
} TLASNode;

float3 transformVector(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V), dot(T->s456, *V), dot(T->s89A, *V));
}

float3 transformPosition(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V) + T->s3, dot(T->s456, *V) + T->s7, dot(T->s89A, *V) + T->sb);
}

void transformRay(Ray *ray, float16 *invTransform)
{
    // do the transform
    ray->D = transformVector(&ray->D, invTransform);
    ray->O = transformPosition(&ray->O, invTransform);
    // update ray direction reciprocals
    ray->rD = (float3)(1.0f / ray->D.x, 1.0f / ray->D.y, 1.0f / ray->D.z);
}

Ray copyRay(Ray *ray)
{
    Ray newRay;
    newRay.O = ray->O;
    newRay.D = ray->D;
    newRay.rD = ray->rD;
    newRay.t = ray->t;
    newRay.barycentric = ray->barycentric;
    newRay.objIdx = ray->objIdx;
    newRay.triIdx = ray->triIdx;
    newRay.traversed = ray->traversed;
    newRay.tested = ray->tested;
    newRay.inside = ray->inside;
    return newRay;
}

void intersectFloor(Ray *ray)
{
    float3 N = (float3)(0, 1, 0); // dist = 1
    float t = -(dot(ray->O, N) + 1) / (dot(ray->D, N));
    if (t < ray->t && t > 0)
        ray->t = t, ray->objIdx = 1;
}

void intersectTri(Ray *ray, const Tri *tri, const int objIdx, uint triIdx)
{
    const float3 vertex0 = (float3)(tri->v0x, tri->v0y, tri->v0z);
    const float3 vertex1 = (float3)(tri->v1x, tri->v1y, tri->v1z);
    const float3 vertex2 = (float3)(tri->v2x, tri->v2y, tri->v2z);
    const float3 edge1 = vertex1 - vertex0;
    const float3 edge2 = vertex2 - vertex0;
    const float3 h = cross(ray->D, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f)
        return; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O - vertex0;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1)
        return;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D, q);
    if (v < 0 || u + v > 1)
        return;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f)
    {
        if (t < ray->t)
            ray->t = min(ray->t, t), ray->objIdx = objIdx, ray->triIdx = triIdx, ray->barycentric = (float2)(u, v);
    }
}

float intersectAABB(const Ray *ray, const float3 bmin, const float3 bmax)
{
    float tx1 = (bmin.x - ray->O.x) * ray->rD.x, tx2 = (bmax.x - ray->O.x) * ray->rD.x;
    float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
    float ty1 = (bmin.y - ray->O.y) * ray->rD.y, ty2 = (bmax.y - ray->O.y) * ray->rD.y;
    tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
    float tz1 = (bmin.z - ray->O.z) * ray->rD.z, tz2 = (bmax.z - ray->O.z) * ray->rD.z;
    tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
    if (tmax >= tmin && tmin < ray->t && tmax > 0)
        return tmin;
    else
        return 1e30f;
}

void intersectBVH(Ray *ray, BVH *bvh, BVHNode *bvhNodes, Tri *triangles, uint *triIdxs, int objIdx)
{
    BVHNode *node = &bvhNodes[bvh->startNodeIdx], *stack[32];
    uint stackPtr = 0;
    while (1)
    {
        if (node->triCount > 0)
        {
            for (uint i = 0; i < node->triCount; i++)
            {
                uint triIdx = triIdxs[node->leftFirst + i];
                Tri *tri = &triangles[triIdx];
                intersectTri(ray, tri, objIdx, triIdx);
            }
            if (stackPtr == 0)
                break;
            else
                node = stack[--stackPtr];
            continue;
        }
        BVHNode *child1 = &bvhNodes[bvh->startNodeIdx + node->leftFirst];
        BVHNode *child2 = &bvhNodes[bvh->startNodeIdx + node->leftFirst + 1];
        float dist1 = intersectAABB(ray, (float3)(child1->aabbMinx, child1->aabbMiny, child1->aabbMinz),
                                    (float3)(child1->aabbMaxx, child1->aabbMaxy, child1->aabbMaxz));
        float dist2 = intersectAABB(ray, (float3)(child2->aabbMinx, child2->aabbMiny, child2->aabbMinz),
                                    (float3)(child2->aabbMaxx, child2->aabbMaxy, child2->aabbMaxz));
        if (dist1 > dist2)
        {
            float d = dist1;
            dist1 = dist2;
            dist2 = d;
            BVHNode *c = child1;
            child1 = child2;
            child2 = c;
        }
        if (dist1 == 1e30f)
        {
            if (stackPtr == 0)
                break;
            else
                node = stack[--stackPtr];
        }
        else
        {
            node = child1;
            if (dist2 != 1e30f)
                stack[stackPtr++] = child2;
        }
    }
}

void intersectBLAS(Ray *ray, int objIdx, float16 invT, BVH *bvh, BVHNode *bvhNodes, Tri *triangles, uint *triIdxs)
{
    Ray tRay = copyRay(ray);
    transformRay(&tRay, &invT);

    // intersectBVH
    intersectBVH(&tRay, bvh, bvhNodes, triangles, triIdxs, objIdx);

    tRay.O = ray->O;
    tRay.D = ray->D;
    tRay.rD = ray->rD;
    ray = &tRay;
}

void intersectTLAS(Ray *ray, TLASNode *tlasNodes, BLAS *blases, BVH *bvhes, BVHNode *bvhNodes, Tri *triangles,
                   uint *triIdxs)
{
    TLASNode *node = &tlasNodes[0], *stack[64];
    uint stackPtr = 0;
    while (1)
    {
        ray->traversed++;
        if (node->leftRight == 0)
        {
            BLAS blas = blases[node->BLAS];
            intersectBLAS(ray, blas.objIdx, blas.invT, &bvhes[blas.bvhIdx], bvhNodes, triangles, triIdxs);
            if (stackPtr == 0)
                break;
            else
                node = stack[--stackPtr];
        }
        TLASNode *child1 = &tlasNodes[node->leftRight & 0xffff];
        TLASNode *child2 = &tlasNodes[node->leftRight >> 16];
        float dist1 = intersectAABB(ray, (float3)(child1->aabbMinx, child1->aabbMiny, child1->aabbMinz),
                                    (float3)(child1->aabbMaxx, child1->aabbMaxy, child1->aabbMaxz));
        float dist2 = intersectAABB(ray, (float3)(child2->aabbMinx, child2->aabbMiny, child2->aabbMinz),
                                    (float3)(child2->aabbMaxx, child2->aabbMaxy, child2->aabbMaxz));
        if (dist1 > dist2)
        {
            float dist = dist1;
            dist1 = dist2;
            dist2 = dist;

            TLASNode *child = child1;
            child1 = child2;
            child2 = child;
        }
        if (dist1 == 1e30f)
        {
            if (stackPtr == 0)
                break;
            else
                node = stack[--stackPtr];
        }
        else
        {
            node = child1;
            if (dist2 != 1e30f)
                stack[stackPtr++] = child2;
        }
    }
}

__kernel void extend(__global Ray *rayBuffer, __global Tri *triBuffer, __global TriEx *triExBuffer,
                     __global uint *triIdxBuffer, __global MeshInstance *meshInstances, __global BVHNode *bvhNodes,
                     __global BVH *bvhes, __global BLAS *blases, __global TLASNode *tlasNodes)
{
    const int index = get_global_id(0);

    Ray ray = rayBuffer[index];

    intersectFloor(&ray);
    intersectTLAS(&ray, tlasNodes, blases, bvhes, bvhNodes, triBuffer, triIdxBuffer);

    rayBuffer[index] = ray;
}