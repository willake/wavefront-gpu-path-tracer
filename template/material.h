#pragma once
#include <memory>

namespace Tmpl8
{
	struct Material
	{
	public:
		Material(const bool isAlbedoOverridden = false)
		{
			bool isLight = false;
			this->albedo = float3(1.0f);
			this->isAlbedoOverridden = isAlbedoOverridden;
			this->reflectivity = 0.0f;
			this->refractivity = 0.0f;
			this->absorption = float3(0);
		}
		//Material(const Material& mat)
		//{
		//	this->type = mat.type;
		//	this->albedo = mat.albedo;
		//	this->isAlbedoOverridden = mat.isAlbedoOverridden;
		//	this->reflectivity = mat.reflectivity;
		//	this->refractivity = mat.refractivity;
		//	this->absorption = float3(0);
		//	textureDiffuse = std::make_unique<Texture>(mat.textureDiffuse.get());
		//}
		float3 GetAlbedo(float2 uv)
		{
			if (textureDiffuse.get() == nullptr)
			{
				return albedo;
			}
			return textureDiffuse->Sample(uv.x, uv.y);
		}
	public:
		bool isLight = false;
		float3 albedo = float3(1.0f);
		bool isAlbedoOverridden = false;
		float reflectivity = 0.0f;
		float refractivity = 0.0f;
		float3 absorption = float3(0.0f);
		std::unique_ptr<Texture> textureDiffuse;
		/*Texture textureMetallic;
		Texture textuteRoughness;*/
	};
}