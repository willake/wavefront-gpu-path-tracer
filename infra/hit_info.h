#pragma once

struct HitInfo
{
public:
	HitInfo() = default;
	HitInfo(float3 n, float2 u, Material* mat) : normal(n), uv(u), material(mat) {}
	float3 normal = float3(0.0f);
	float2 uv = float2(0.0f);
	Material* material;
};