typedef struct __attribute__((aligned(128)))
{
    float3 O, D, rD, E; // 64 bytes
    float t;            // 4 bytes
    int pixelIdx;       // 4 bytes
} ShadowRay;            // total 72 bytes

ShadowRay copyRay(ShadowRay *ray)
{
    ShadowRay newRay;
    newRay.O = ray->O;
    newRay.D = ray->D;
    newRay.rD = ray->rD;
    newRay.E = ray->E;
    newRay.t = ray->t;
    newRay.pixelIdx = ray->pixelIdx;
    return newRay;
}

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
    float uv0x, uv0y;
    float uv1x, uv1y;
    float uv2x, uv2y;
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
} TLASNode; // 32 bytes in total

typedef struct
{
    float16 T;                    // 64 bytes
    float16 invT;                 // 64 bytes
    float colorx, colory, colorz; // 12 bytes
    float size;                   // 4 bytes
    int objIdx;                   // 4 bytes
} Light;                          // 148 bytes in total

float3 transformVector(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V), dot(T->s456, *V), dot(T->s89A, *V));
}

float3 transformPosition(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V) + T->s3, dot(T->s456, *V) + T->s7, dot(T->s89A, *V) + T->sb);
}

void transformRay(ShadowRay *ray, float16 *invTransform)
{
    // do the transform
    ray->D = transformVector(&ray->D, invTransform);
    ray->O = transformPosition(&ray->O, invTransform);
    // update ray direction reciprocals
    ray->rD = (float3)(1.0f / ray->D.x, 1.0f / ray->D.y, 1.0f / ray->D.z);
}

bool isOccludedLight(const ShadowRay *ray, Light *light)
{
    float16 invT = light->invT;
    float size = light->size * 0.5f;
    const float Oy = invT[4] * ray->O.x + invT[5] * ray->O.y + invT[6] * ray->O.z + invT[7];
    const float Dy = invT[4] * ray->D.x + invT[5] * ray->D.y + invT[6] * ray->D.z;
    const float t = Oy / -Dy;
    if (t < ray->t && t > 0)
    {
        const float Ox = invT[0] * ray->O.x + invT[1] * ray->O.y + invT[2] * ray->O.z + invT[3];
        const float Oz = invT[8] * ray->O.x + invT[9] * ray->O.y + invT[10] * ray->O.z + invT[11];
        const float Dx = invT[0] * ray->D.x + invT[1] * ray->D.y + invT[2] * ray->D.z;
        const float Dz = invT[8] * ray->D.x + invT[9] * ray->D.y + invT[10] * ray->D.z;
        const float Ix = Ox + t * Dx, Iz = Oz + t * Dz;
        if (Ix > -size && Ix < size && Iz > -size && Iz < size)
        {
            return true;
        }
    }
    return false;
}

bool isOccludedFloor(const ShadowRay *ray)
{
    float3 N = (float3)(0, 1, 0); // dist = 1
    float t = -(dot(ray->O, N) + 1) / (dot(ray->D, N));
    if (t < ray->t && t > 0)
    {
        return true;
    }
    return false;
}

bool isOccludedTri(const ShadowRay *ray, const Tri *tri)
{
    const float3 vertex0 = (float3)(tri->v0x, tri->v0y, tri->v0z);
    const float3 vertex1 = (float3)(tri->v1x, tri->v1y, tri->v1z);
    const float3 vertex2 = (float3)(tri->v2x, tri->v2y, tri->v2z);
    const float3 edge1 = vertex1 - vertex0;
    const float3 edge2 = vertex2 - vertex0;
    const float3 h = cross(ray->D, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f)
        return false; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O - vertex0;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1)
        return false;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D, q);
    if (v < 0 || u + v > 1)
        return false;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f)
    {
        if (t < ray->t)
        {
            return true;
        }
    }
    return false;
}

float intersectAABB(const ShadowRay *ray, const float3 bmin, const float3 bmax)
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

bool isOccludedBVH(const ShadowRay *ray, BVH *bvh, MeshInstance *meshIns, BVHNode *bvhNodes, Tri *triangles,
                   uint *triIdxs)
{
    BVHNode *node = &bvhNodes[bvh->startNodeIdx], *stack[32];
    uint stackPtr = 0;
    while (1)
    {
        if (node->triCount > 0)
        {
            for (uint i = 0; i < node->triCount; i++)
            {
                uint triIdx = triIdxs[meshIns->triStartIdx + node->leftFirst + i];
                Tri *tri = &triangles[triIdx];
                if (isOccludedTri(ray, tri))
                {
                    return true;
                }
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
    return false;
}

bool isOccludedBLAS(ShadowRay *ray, BVH *bvh, MeshInstance *meshIns, BVHNode *bvhNodes, Tri *triangles, uint *triIdxs,
                    float16 invT)
{
    ShadowRay tRay = copyRay(ray);
    transformRay(&tRay, &invT);
    // intersectBVH
    return isOccludedBVH(&tRay, bvh, meshIns, bvhNodes, triangles, triIdxs);
}

bool isOccludedTLAS(ShadowRay *ray, TLASNode *tlasNodes, BLAS *blases, BVH *bvhes, MeshInstance *meshInstances,
                    BVHNode *bvhNodes, Tri *triangles, uint *triIdxs)
{
    TLASNode *node = &tlasNodes[0], *stack[32];
    uint stackPtr = 0;
    while (1)
    {
        if (node->leftRight == 0) // isLeaf
        {
            // ray->traversed += 100;
            BLAS blas = blases[node->BLAS];
            if (isOccludedBLAS(ray, &bvhes[blas.bvhIdx], &meshInstances[blas.bvhIdx], bvhNodes, triangles, triIdxs,
                               blas.invT))
            {
                return true;
            }
            if (stackPtr == 0)
                break;
            else
                node = stack[--stackPtr];
            continue;
        }
        TLASNode *child1 = &tlasNodes[node->leftRight & 0xffff];
        TLASNode *child2 = &tlasNodes[node->leftRight >> 16];
        float dist1 = intersectAABB(ray, (float3)(child1->aabbMinx, child1->aabbMiny, child1->aabbMinz),
                                    (float3)(child1->aabbMaxx, child1->aabbMaxy, child1->aabbMaxz));
        float dist2 = intersectAABB(ray, (float3)(child2->aabbMinx, child2->aabbMiny, child2->aabbMinz),
                                    (float3)(child2->aabbMaxx, child2->aabbMaxy, child2->aabbMaxz));
        if (dist1 > dist2)
        {
            float d = dist1;
            dist1 = dist2;
            dist2 = d;

            TLASNode *c = child1;
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
    return false;
}

__kernel void connect(__global float4 *Ts, __global float4 *Es, __global ShadowRay *rayBuffer, __global Tri *triBuffer,
                      __global uint *triIdxBuffer, __global BVHNode *bvhNodes, __global BVH *bvhes,
                      __global BLAS *blases, __global TLASNode *tlasNodes, __global MeshInstance *meshInstances,
                      __global Light *lights, uint lightCount)
{
    const int index = get_global_id(0);

    ShadowRay ray = rayBuffer[index];

    for (int i = 0; i < lightCount; i++)
    {
        if (isOccludedLight(&ray, &lights[i]))
        {
            return;
        }
    }

    if (isOccludedFloor(&ray))
    {
        return;
    }

    if (isOccludedTLAS(&ray, tlasNodes, blases, bvhes, meshInstances, bvhNodes, triBuffer, triIdxBuffer))
    {
        return;
    }

    Es[ray.pixelIdx] += Ts[ray.pixelIdx] * (float4)(ray.E.x, ray.E.y, ray.E.z, 0);
}