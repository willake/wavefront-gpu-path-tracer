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
    uint objIdx, matIdx, bvhIdx;
    float16 T, invT;
    float aabbMinx, aabbMiny, aabbMinz;
    float aabbMaxx, aabbMaxy, aabbMaxz;
} BLAS;

typedef struct
{
    float3 normal;
    float2 uv;
    int matIdx;
} HitInfo;

float3 transformVector(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V), dot(T->s456, *V), dot(T->s89A, *V));
}

float3 getBLASNormal(TriEx *triExs, const int triIdx, const float2 barycentric, float16 T)
{
    TriEx *ex = &triExs[triIdx];
    float3 n0 = (float3)(ex->n0x, ex->n0y, ex->n0z);
    float3 n1 = (float3)(ex->n1x, ex->n1y, ex->n1z);
    float3 n2 = (float3)(ex->n2x, ex->n2y, ex->n2z);
    float3 N = (1 - barycentric.x - barycentric.y) * n0 + barycentric.x * n1 + barycentric.y * n2;
    transformVector(&N, &T);
    return normalize(N);
}

float2 getBLASUV(TriEx *triExs, const int triIdx, const float2 barycentric)
{
    TriEx *ex = &triExs[triIdx];
    float2 uv0 = (float2)(ex->uv0x, ex->uv0y);
    float2 uv1 = (float2)(ex->uv1x, ex->uv1y);
    float2 uv2 = (float2)(ex->uv2x, ex->uv2y);
    return (1 - barycentric.x - barycentric.y) * uv0 + barycentric.x * uv1 + barycentric.y * uv2;
}

float3 getFloorNormal()
{
    return (float3)(0, 1, 0);
}
float2 getFloorUV(const float3 I)
{
    float u = I.x;
    float v = I.z;

    const float invto = 1.0f / 5.12f;

    u *= invto;
    v *= invto;

    u = u - floor(u);
    v = v - floor(v);

    // Return the UV coordinates
    return (float2)(u, v);
}

uint getSkyColor(Ray *ray, uint *pixels, uint width, uint height)
{
    if (!pixels)
        return 0;

    float phi = atan2(-ray->D.z, ray->D.x);
    uint u = (uint)(width * (phi > 0 ? phi : (phi + 2 * M_PI_F)) * M_2_PI_F - 0.5f);
    uint v = (uint)(height * acos(ray->D.y) * M_1_PI_F - 0.5f);
    uint skyIdx = (u + v * width) % (width * height);

    return pixels[skyIdx];
}

uint sample(uint *pixels, float2 uv, uint width, uint height)
{

    float u = clamp(uv.x, 0.0f, 1.0f);
    float v = 1 - clamp(uv.y, 0.0f, 1.0f);

    uint x = (uint)(u * width);
    uint y = (uint)(v * height);

    x = clamp(x, (uint)0, width - 1);
    y = clamp(y, (uint)0, height - 1);

    uint index = x + y * width;

    return pixels[index];
}

HitInfo getHitInfo(const Ray *ray, TriEx *triExs, BLAS *blases, const float3 I)
{
    HitInfo hitInfo;
    if (ray->objIdx == 1)
    {
        hitInfo.normal = getFloorNormal();
        hitInfo.uv = getFloorUV(I);
    }
    else
    {
        BLAS *blas = &blases[ray->objIdx - 2];
        hitInfo.normal = getBLASNormal(triExs, ray->triIdx, ray->barycentric, blas->T);
        hitInfo.uv = getBLASUV(triExs, ray->triIdx, ray->barycentric);
        hitInfo.matIdx = blas->matIdx;
    }

    return hitInfo;
}

__kernel void shade(__global uint *accumulator, __global Ray *rayBuffer, __global uint *skydomePixels,
                    uint skydomeWidth, uint skydomeHeight, __global uint *floorPixels, __global TriEx *triExs,
                    __global BLAS *blases)
{
    // get ray id
    const int index = get_global_id(0);

    Ray ray = rayBuffer[index];

    if (ray.objIdx == -1)
    {
        accumulator[index] = getSkyColor(&ray, skydomePixels, skydomeWidth, skydomeHeight);
        return;
    }

    float3 I = ray.O + ray.D * ray.t;
    HitInfo hitInfo = getHitInfo(&ray, triExs, blases, I);
    float3 N = hitInfo.normal;
    float2 uv = hitInfo.uv;

    if (ray.objIdx == 1)
    {
        accumulator[index] = sample(floorPixels, hitInfo.uv, 512, 512);
    }

    if (ray.traversed > 2)
    {
        accumulator[index] = 9527;
    }

    if (ray.traversed > 4)
    {
        accumulator[index] = 24601;
    }

    if (ray.traversed > 5)
    {
        accumulator[index] = 24601496;
    }
}