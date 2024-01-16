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
    bool lastSpecular;  // 1 bytes
} Ray;                  // total 82 bytes

inline float4 getTraverseCountColor(int traversed, int peak)
{
    const float invMax = 1 / 255.f;
    const float4 green = (float4)(179 * invMax, 255 * invMax, 174 * invMax, 0.0f);
    const float4 red = (float4)(255 * invMax, 50 * invMax, 50 * invMax, 0.0f);

    if (peak < 10)
        return green;

    traversed = clamp(traversed, 0, peak);
    float blend = traversed / (float)peak;

    float r = green.x + blend * (red.x - green.x);
    float g = green.y + blend * (red.y - green.y);
    float b = green.z + blend * (red.z - green.z);

    return (float4)(r, g, b, 0);
}

__kernel void renderTraversal(__global float4 *Es, __global Ray *rayBuffer)
{
    const int index = get_global_id(0);

    Ray ray = rayBuffer[index];

    Es[index] = getTraverseCountColor(ray.traversed, 100);
}