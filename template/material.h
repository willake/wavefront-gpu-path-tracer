#pragma once
#include <memory>

namespace Tmpl8
{
__declspec(align(32)) struct GPUMaterial
{
    GPUMaterial(){};
    float3 albedo = float3(1.0f);     // 12 bytes
    float3 absorption = float3(0.0f); // 12 bytes
    float reflectivity = 0.0f;        // 4 bytes
    float refractivity = 0.0f;        // 4 bytes
                                      // 32 bytes in total
};
struct Material
{
  public:
    Material(const bool albedoOverridden = false)
    {
        isLight = false;
        albedo = float3(1.0f);
        isAlbedoOverridden = albedoOverridden;
        absorption = float3(0);
        roughness = 0.6f;
        metalness = 0.5f;
        transmittance = 0;
    }

    float3 GetAlbedo(float2 uv)
    {
        /*if (textureDiffuse.get() == nullptr)
        {
            return albedo;
        }*/
        return textureDiffuse.Sample(uv.x, uv.y);
    }

    float3 Evaluate(float3 L, float3 N, float3 V, float3 baseColor)
    {
        float3 H = normalize(V + L); // H = normalize(V + L);
        float NdotH = dot(N, H), NdotV = dot(N, V), NdotL = dot(N, L), VdotH = dot(V, H);

        float3 f0 = float3(0.04f);
        f0 = mix(f0, baseColor, metalness);

        float3 F = FresnelSchlick(VdotH, f0);
        float D = DistributionGGX(NdotH, roughness);
        float G = GeometrySmith(NdotV, NdotL, roughness);

        float3 spec = F * G * D / 4.0f * max(0.00001f, NdotL) * max(0.00001f, NdotV);

        baseColor *= float3(1.0) - F;

        baseColor *= (1.0 - metalness);

        float3 diffuse = baseColor * INVPI * NdotL;

        return diffuse + spec;
    }

  public:
    bool isLight = false;
    float3 albedo = float3(1.0f);
    bool isAlbedoOverridden = false;
    float roughness;
    float metalness;
    float transmittance;
    float3 absorption = float3(0.0f);
    Texture textureDiffuse;
    float reflectivity = 0.0f;
    float refractivity = 0.0f;
    // std::unique_ptr<Texture> textureDiffuse;
    /*Texture textureMetallic;
    Texture textuteRoughness;*/
};
} // namespace Tmpl8