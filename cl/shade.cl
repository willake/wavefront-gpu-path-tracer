#define EPSILON 0.001f
// random numbers: seed using WangHash((threadidx+1)*17), then use RandomInt / RandomFloat
uint WangHash(uint s)
{
    s = (s ^ 61) ^ (s >> 16), s *= 9, s = s ^ (s >> 4), s *= 0x27d4eb2d, s = s ^ (s >> 15);
    return s;
}
uint RandomInt(uint *s)
{
    *s ^= *s << 13, *s ^= *s >> 17, *s ^= *s << 5;
    return *s;
}
float RandomFloat(uint *s)
{
    return RandomInt(s) * 2.3283064365387e-10f; /* = 1 / (2^32-1) */
}

// Define ray
// __attribute__((aligned(64)))
typedef struct __attribute__((aligned(128)))
{
    float3 O, D, rD;    // 48 bytes
    float2 barycentric; // 8 bytes
    float t;            // 4 bytes
    int objIdx;         // 4 bytes
    int triIdx;         // 4 bytes
    int pixelIdx;       // 4 bytes
    int traversed;      // 4 bytes
    int tested;         // 4 bytes
    bool inside;        // 1 bytes
} Ray;                  // total 81 bytes

Ray GenerateRay(const float3 origin, const float3 direction, const int pixelIdx)
{
    Ray ray;
    ray.O = origin, ray.D = direction, ray.t = 1e34f;
    ray.rD = (float3)(1 / ray.D.x, 1 / ray.D.y, 1 / ray.D.z);
    ray.objIdx = -1;
    ray.pixelIdx = pixelIdx;
    return ray;
}

typedef struct __attribute__((aligned(128)))
{
    float3 O, D, rD, E; // 64 bytes
    float t;            // 4 bytes
    int pixelIdx;       // 4 bytes
} ShadowRay;            // total 72 bytes

ShadowRay GenerateShadowRay(const float3 origin, const float3 direction, const float distance, const int pixelIdx)
{
    ShadowRay shadowRay;
    shadowRay.O = origin, shadowRay.D = direction, shadowRay.t = distance;
    // ray.rD = (float3)(1 / ray.D.x, 1 / ray.D.y, 1 / ray.D.z);
    // ray.pixelIdx = pixelIdx;
    return shadowRay;
}

typedef struct __attribute__((aligned(64)))
{
    float2 uv0, uv1, uv2; // 24 bytes
    float n0x, n0y, n0z;  // 12 bytes
    float n1x, n1y, n1z;  // 12 bytes
    float n2x, n2y, n2z;  // 12 bytes
    float dummy;          // 4 bytes
} TriEx;                  // total 64 bytes

typedef struct
{
    uint objIdx, matIdx, bvhIdx;        // 12 bytes
    float16 T, invT;                    // 128 bytes
    float aabbMinx, aabbMiny, aabbMinz; // 12 bytes
    float aabbMaxx, aabbMaxy, aabbMaxz; // 12 bytes
} BLAS;

typedef struct __attribute__((aligned(32)))
{
    float albedox, albedoy, albedoz;             // 12 bytes
    float absorptionx, absorptiony, absorptionz; // 12 bytes
    float reflectivity;                          // 4 bytes
    float refractivity;                          // 4 bytes
} Material;                                      // 32 bytes

typedef struct __attribute__((aligned(16)))
{
    uint width;
    uint height;
    uint startIdx;
    uint dummy;
} Texture;

typedef struct
{
    float16 T;                    // 64 bytes
    float16 invT;                 // 64 bytes
    float colorx, colory, colorz; // 12 bytes
    float size;                   // 4 bytes
    int objIdx;                   // 4 bytes
} Light;                          // 148 bytes in total

typedef struct
{
    float3 normal;
    float2 uv;
    int matIdx;
} HitInfo;

float3 diffusereflection(float3 *N, uint *seed)
{
    float3 R;
    do
    {
        R = (float3)(RandomFloat(seed) * 2 - 1, RandomFloat(seed) * 2 - 1, RandomFloat(seed) * 2 - 1);
    } while (dot(R, R) > 1);
    if (dot(R, *N) < 0)
        R *= -1.0f;
    return normalize(R);
}

uint RGB32FtoRGB8(float3 c)
{
    int r = (int)(min(c.x, 1.f) * 255);
    int g = (int)(min(c.y, 1.f) * 255);
    int b = (int)(min(c.z, 1.f) * 255);
    return (r << 16) + (g << 8) + b;
}

float3 RGB8toRGB32F(uint c)
{
    float s = 1 / 256.0f;
    int r = (c >> 16) & 255;
    int g = (c >> 8) & 255;
    int b = c & 255;
    return (float3)(r * s, g * s, b * s);
}

float3 transformVector(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V), dot(T->s456, *V), dot(T->s89A, *V));
}

float3 transformPosition(float3 *V, float16 *T)
{
    return (float3)(dot(T->s012, *V) + T->s3, dot(T->s456, *V) + T->s7, dot(T->s89A, *V) + T->sb);
}

float3 getLightNormal(Light *light)
{
    float3 N = (float3)(0, -1, 0);
    return transformVector(&N, &light->T);
}

float3 getBLASNormal(TriEx *triExs, const int triIdx, const float2 barycentric, float16 T)
{
    TriEx *ex = &triExs[triIdx];
    float3 n0 = (float3)(ex->n0x, ex->n0y, ex->n0z);
    float3 n1 = (float3)(ex->n1x, ex->n1y, ex->n1z);
    float3 n2 = (float3)(ex->n2x, ex->n2y, ex->n2z);
    float3 N = (1 - barycentric.x - barycentric.y) * n0 + barycentric.x * n1 + barycentric.y * n2;
    N = transformVector(&N, &T);
    return normalize(N);
}

float2 getBLASUV(TriEx *triExs, const int triIdx, const float2 barycentric)
{
    TriEx *ex = &triExs[triIdx];
    float2 uv0 = ex->uv0;
    float2 uv1 = ex->uv1;
    float2 uv2 = ex->uv2;
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

/* For Debugging */
float3 getEdgeColor(float2 barycentric)
{
    if (fabs(barycentric.x) < 0.03f || fabs(barycentric.x - 1) < 0.03f || fabs(barycentric.y) < 0.03f ||
        fabs(barycentric.y - 1) < 0.03f)
    {
        return (float3)(0, 0, 0);
    }
    else
    {
        return (float3)(1);
    }
}

float3 getSkyColor(Ray *ray, uint *pixels, uint width, uint height)
{
    if (!pixels)
        return 0;

    float phi = atan2(-ray->D.z, ray->D.x);
    uint u = (uint)(width * (phi > 0 ? phi : (phi + 2 * M_PI_F)) * M_2_PI_F - 0.5f);
    uint v = (uint)(height * acos(ray->D.y) * M_1_PI_F - 0.5f);
    uint skyIdx = (u + v * width) % (width * height);

    return RGB8toRGB32F(pixels[skyIdx]);
}

float3 sample(uint *pixels, uint startIdx, float2 uv, uint width, uint height)
{

    float u = clamp(uv.x, 0.0f, 1.0f);
    float v = 1 - clamp(uv.y, 0.0f, 1.0f);

    uint x = (uint)(u * width);
    uint y = (uint)(v * height);

    x = clamp(x, (uint)0, width - 1);
    y = clamp(y, (uint)0, height - 1);

    uint index = x + y * width;

    return RGB8toRGB32F(pixels[startIdx + index]);
}

HitInfo getHitInfo(const Ray *ray, TriEx *triExs, BLAS *blases, Light *lights, const float3 I)
{
    HitInfo hitInfo;
    if (ray->objIdx >= 900)
    {
        hitInfo.normal = getLightNormal(&lights[ray->objIdx - 900]);
        hitInfo.uv = (float2)(0, 0);
    }
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

    if (ray->objIdx < 900 && dot(hitInfo.normal, ray->D) > 0)
        hitInfo.normal = -hitInfo.normal;

    return hitInfo;
}

float3 getAlbedo(__global uint *floorPixels, __global BLAS *blases, __global uint *texturePixels,
                 __global Texture *textures, __global Light *lights, int objIdx, float2 uv)
{
    if (objIdx >= 900)
    {
        Light *light = &lights[objIdx - 900];
        return (float3)(light->colorx, light->colory, light->colorz);
    }
    else if (objIdx == 1)
    {
        return sample(floorPixels, 0, uv, 512, 512);
    }
    else
    {
        BLAS *blas = &blases[objIdx - 2];
        Texture *texture = &textures[blas->matIdx];
        return sample(texturePixels, texture->startIdx, uv, texture->width, texture->height);
    }
}

float3 randomPointOnLight(Light *lights, uint lightCount, uint *seed, uint *lightIdx)
{
    float r0 = RandomFloat(seed);
    float r1 = RandomFloat(seed);
    *lightIdx = (uint)(r0 * lightCount);
    Light *light = &lights[*lightIdx];
    // renormalize r0 for reuse
    float stratum = *lightIdx * 0.25f;
    float r2 = (r0 - stratum) / (1 - stratum);
    // get a random position on the select quad
    const float size = light->size * 0.5f;
    float3 corner1 = transformPosition(&(float3)(-size, 0, -size), &light->T);
    float3 corner2 = transformPosition(&(float3)(size, 0, -size), &light->T);
    float3 corner3 = transformPosition(&(float3)(-size, 0, size), &light->T);
    return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
}

ShadowRay directionIllumination(Light *lights, uint lightCount, uint *seed, float3 I, float3 N, float3 brdf,
                                int pixelIdx)
{
    uint lightIdx;
    float3 randomLightPos = randomPointOnLight(lights, lightCount, seed, &lightIdx);
    Light *light = &lights[lightIdx];
    // lighting related information
    float3 L = randomLightPos - I;
    float dist = length(L);
    L = normalize(L);
    float ndotl = dot(N, L);
    float nldotl = dot(getLightNormal(light), -L);
    float A = light->size * light->size;

    ShadowRay shadowRay = GenerateShadowRay(I + L * EPSILON, L, dist - 2 * EPSILON, pixelIdx);
    if (ndotl > 0 && nldotl > 0)
    {
        float solidAngle = (nldotl * A) / (dist * dist);

        // ray.E = (float3)(light->colorx, light->colory, light->colorz) * solidAngle * brdf * ndotl * lightCount;
    }
    return shadowRay;
}

__kernel void shade(__global float4 *pixels, __global Ray *rayBuffer, __global uint *seeds,
                    __global uint *skydomePixels, uint skydomeWidth, uint skydomeHeight, __global uint *floorPixels,
                    __global TriEx *triExs, __global BLAS *blases, __global Material *materials,
                    __global uint *texturePixels, __global Texture *textures, __global Light *lights, uint lightCount,
                    __global Ray *extensionrayBuffer, __global ShadowRay *shadowrayBuffer,
                    __global uint *extensionrayCounter, __global uint *shadowrayCounter, uint depth)
{
    // get ray id
    const int index = get_global_id(0);

    Ray ray = rayBuffer[index];
    const int pixelIdx = ray.pixelIdx;
    uint seed = seeds[pixelIdx];

    if (ray.objIdx == -1)
    {
        float3 skyColor = getSkyColor(&ray, skydomePixels, skydomeWidth, skydomeHeight);
        pixels[pixelIdx] *= (float4)(skyColor.x, skyColor.y, skyColor.z, 1);
        return;
    }

    float3 I = ray.O + ray.D * ray.t;
    HitInfo hitInfo = getHitInfo(&ray, triExs, blases, lights, I);
    float3 N = hitInfo.normal;
    float2 uv = hitInfo.uv;
    float3 albedo = getAlbedo(floorPixels, blases, texturePixels, textures, lights, ray.objIdx, uv);
    float3 brdf = albedo * M_1_PI_F;
    /* visualize triangle */
    // accumulator[index] = ray.triIdx * 10;
    // return;
    /* visualize normal */                   // float3 color = (N + 1) * 0.5f;
    /* visualize uv */                       // float3 color = (float3)(uv.x, uv.y, 0);
    /* visualize visualize triangle edges */ // float3 color = getEdgeColor(ray.barycentric);
    /* debug */                              // accumulator[index] = RGB32FtoRGB8(color); return;

    // objIdx >= 900 is a light
    if (ray.objIdx >= 900)
    {
        pixels[pixelIdx] = (float4)(0);
        return;
    }

    BLAS blas = blases[ray.objIdx - 2];
    Material material = materials[blas.matIdx];

    float3 out_radiance = (float3)(0);
    float reflectivity = material.reflectivity;
    float refractivity = material.reflectivity;

    float3 medium_scale = (float3)(1);
    if (ray.inside)
    {
        float3 absorption = (float3)(material.absorptionx, material.absorptiony, material.absorptionz);
        medium_scale = exp(absorption * -ray.t);
    }

    float r = RandomFloat(&seed);
    if (r < reflectivity) // handle pure speculars
    {
        // generate extend ray
    }
    else if (r < reflectivity + refractivity) // handle dielectrics
    {
        // generate extend ray
    }
    else // diffuse surface
    {
        float3 R = diffusereflection(&N, &seed);
        // compute diffuse
        pixels[pixelIdx] *= (float4)(medium_scale.x, medium_scale.y, medium_scale.z, 0) *
                            (float4)(brdf.x, brdf.y, brdf.z, 0) * 2 * M_PI_F * dot(R, N);
        // generate extend ray
        uint ei = atomic_inc(extensionrayCounter);
        extensionrayBuffer[ei] = GenerateRay(I + R * EPSILON, R, pixelIdx);

        uint si = atomic_inc(shadowrayCounter);
        shadowrayBuffer[si] = directionIllumination(lights, lightCount, &seed, I, N, brdf, pixelIdx);
        // if (depth == 0)
        //     pixels[pixelIdx] *= (float4)(brdf.x, brdf.y, brdf.z, 0);
    }
    // pixels[pixelIdx] *= (float4)(albedo.x, albedo.y, albedo.z, 1);
}