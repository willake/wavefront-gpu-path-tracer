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
        reflectivity = 0.0f;
        refractivity = 0.0f;
        absorption = float3(0);
    }

    float3 GetAlbedo(float2 uv)
    {
        /*if (textureDiffuse.get() == nullptr)
        {
            return albedo;
        }*/
        return textureDiffuse.Sample(uv.x, uv.y);
    }

  public:
    bool isLight = false;
    float3 albedo = float3(1.0f);
    bool isAlbedoOverridden = false;
    float reflectivity = 0.0f;
    float refractivity = 0.0f;
    float3 absorption = float3(0.0f);
    Texture textureDiffuse;
    // std::unique_ptr<Texture> textureDiffuse;
    /*Texture textureMetallic;
    Texture textuteRoughness;*/
};
} // namespace Tmpl8