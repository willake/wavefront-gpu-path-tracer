#pragma once

#include <string>
#include <unordered_map>

struct Light
{
    mat4 T;       // 64 bytes
    mat4 invT;    // 64 bytes
    float3 color; // 12 bytes
    float size;   // 4 bytes
    int objIdx;   // 4 bytes
    float3 GetNormal(const float3 I) const
    {
        return TransformVector(float3(0, -1, 0), T);
        // return float3(-T.cell[1], -T.cell[5], -T.cell[9]);
    }
}; // 152 bytes in total

struct Tri
{
    Tri()
    {
        vertex0 = float3(0), vertex1 = float3(0), vertex2 = float3(0);
        centroid = float3(0);
    };
    Tri(float3 v0, float3 v1, float3 v2, float3 c) : vertex0(v0), vertex1(v1), vertex2(v2), centroid(c){};
    float3 vertex0, vertex1, vertex2; // 36 bytes
    float3 centroid;                  // 12 bytes

    aabb GetBounds()
    {
        aabb bounds;
        bounds.Grow(vertex0);
        bounds.Grow(vertex1);
        bounds.Grow(vertex2);
        return bounds;
    }
};

__declspec(align(64)) struct TriEx
{
    TriEx()
    {
        normal0 = float3(0), normal1 = float3(0), normal2 = float3(0);
        uv0 = float2(0), uv1 = float2(0), uv2 = float2(0);
        dummy = 0;
    };
    TriEx(float3 n0, float3 n1, float3 n2, float2 u0, float2 u1, float2 u2)
        : normal0(n0), normal1(n1), normal2(n2), uv0(u0), uv1(u1), uv2(u2){};
    float2 uv0, uv1, uv2;             // 24 bytes
    float3 normal0, normal1, normal2; // 36 bytes
    float dummy;                      // 4 bytes 64 in total
};

struct Vertex
{
    float3 position{};
    float3 normal{};
    float2 uv{};

    bool operator==(const Vertex &other) const
    {
        return position == other.position && normal == other.normal && uv == other.uv;
    }
};

namespace std
{
template <class T> inline void hash_combine(std::size_t &s, const T &v)
{
    std::hash<T> h;
    s ^= h(v) + 0x9e3779b9 + (s << 6) + (s >> 2);
}

template <> struct hash<float2>
{
    size_t operator()(float2 const &v) const
    {
        std::size_t res = 0;
        hash_combine(res, v.x);
        hash_combine(res, v.y);
        return res;
    }
};

template <> struct hash<float3>
{
    size_t operator()(float3 const &v) const
    {
        std::size_t res = 0;
        hash_combine(res, v.x);
        hash_combine(res, v.y);
        hash_combine(res, v.z);
        return res;
    }
};

template <> struct hash<Vertex>
{
    size_t operator()(Vertex const &vertex) const
    {
        return ((hash<float3>()(vertex.position) ^ (hash<float3>()(vertex.normal) << 1)) >> 1) ^
               (hash<float2>()(vertex.uv) << 1);
    }
};
} // namespace std

namespace Tmpl8
{

inline float3 GenerateRandomUnitVector()
{
    // Generate random angles
    const float theta = acosf(2 * RandomFloat() - 1) - (PI / 2);
    const float phi = 2 * PI * RandomFloat();

    // Calculate the corresponding unit vector components
    const float x = std::sin(phi) * std::cos(theta);
    const float y = std::sin(phi) * std::sin(theta);
    const float z = std::cos(phi);

    return float3(x, y, z);
}

inline float3 GetTraverseCountColor(int traversed, int peak)
{
    const float invMax = 1 / 255.f;
    const float3 green(179 * invMax, 255 * invMax, 174 * invMax);
    const float3 red(255 * invMax, 50 * invMax, 50 * invMax);

    if (peak < 10)
        return green;

    traversed = clamp(traversed, 0, peak);
    float blend = traversed / (float)peak;

    float r = green.x + blend * (red.x - green.x);
    float g = green.y + blend * (red.y - green.y);
    float b = green.z + blend * (red.z - green.z);

    return float3(r, g, b);
}

inline float3 GetDepthColor(int current, int max)
{
    current = clamp(current, 0, max);
    float blend = current / (float)max;

    return float3(blend, 1 - blend, 0);
}

inline void ImGuiFloat3(const std::string label, float3 &value)
{
    ImGui::Text("%s", label.c_str());
    ImGui::PushID(label.c_str());
    ImGui::SameLine(); // Move cursor to the same line as the label
    ImGui::SetNextItemWidth(50.0f);
    ImGui::InputFloat("X", &value.x, 0.0f, 0.0f, "%.2f");
    ImGui::SetNextItemWidth(50.0f);
    ImGui::SameLine();
    ImGui::InputFloat("Y", &value.y, 0.0f, 0.0f, "%.2f");
    ImGui::SetNextItemWidth(50.0f);
    ImGui::SameLine();
    ImGui::InputFloat("Z", &value.z, 0.0f, 0.0f, "%.2f");
    ImGui::PopID();
}

uint64_t NowInMicro()
{
    return std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::high_resolution_clock::now().time_since_epoch())
        .count();
}

float3 RGB8toRGB32F(uint c)
{
    float rgbScale = 1 / 256.0f;
    float r = ((c >> 16) & 0xFF) * rgbScale;
    float g = ((c >> 8) & 0xFF) * rgbScale;
    float b = (c & 0xFF) * rgbScale;
    return float3(r, g, b);
}

const float InvGamma = 1.0f / 2.2f;

float4 GammaCorrection(float4 v)
{
    return float4(powf(v.x, InvGamma), powf(v.y, InvGamma), powf(v.z, InvGamma), 0);
}
const float Deg2Red = (PI * 2) / 360.0f;

float3 FresnelSchlick(float cosTheta, float3 f0)
{
    return f0 + (1.0 - f0) * pow(1.0f - cosTheta, 5);
}

float DistributionGGX(float NdotH, float roughness)
{
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float NdotH2 = NdotH * NdotH;
    float b = (NdotH2 * (alpha2 - 1.0f) + 1.0f);
    return alpha2 * INVPI / (b * b);
}

float G1_GGX_Schlick(float NdotV, float roughness)
{
    float alpha = roughness * roughness;
    float k = alpha / 2.0f;
    return max(NdotV, 0.001f) / (NdotV * (1.0f - k) + k);
}

float GeometrySmith(float NdotV, float NdotL, float roughness)
{
    return G1_GGX_Schlick(NdotL, roughness) * G1_GGX_Schlick(NdotV, roughness);
}

float SurvivalProb(float3 color)
{
    return clamp(max(color.x, max(color.y, color.z)), 0.1, 0.9);
}

float3 mix(float3 x, float3 y, float weight)
{
    return x * (1.0f - weight) + y * weight;
}
} // namespace Tmpl8