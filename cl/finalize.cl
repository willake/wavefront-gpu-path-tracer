#define INVGAMMA 1.0f / 2.2f
uint RGB32FtoRGB8(float4 c)
{
    int r = (int)(min(c.x, 1.f) * 255);
    int g = (int)(min(c.y, 1.f) * 255);
    int b = (int)(min(c.z, 1.f) * 255);
    return (r << 16) + (g << 8) + b;
}

float4 RGB8toRGB32F(uint c)
{
    float s = 1 / 256.0f;
    int r = (c >> 16) & 255;
    int g = (c >> 8) & 255;
    int b = c & 255;
    return (float4)(r * s, g * s, b * s, 0);
}

float4 gammaCorrection(float4 v)
{
    return (float4)(pow(v.x, INVGAMMA), pow(v.y, INVGAMMA), pow(v.z, INVGAMMA), 0);
}

__kernel void finalize(__global uint *finalPixels, __global float4 *accumulator, __global float4 *pixels, float scale,
                       uint isDebugging)
{
    const int index = get_global_id(0);

    if (isDebugging)
    {
        finalPixels[index] = RGB32FtoRGB8(pixels[index]);
        return;
    }

    accumulator[index] += pixels[index];
    float4 finalColor = accumulator[index] * scale;
    finalColor = gammaCorrection(finalColor);
    finalPixels[index] = RGB32FtoRGB8(finalColor);
}