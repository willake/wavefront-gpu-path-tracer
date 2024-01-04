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

struct GPUHitInfo
{
	GPUHitInfo() = default;
	GPUHitInfo(float3 n, float2 u, int idx) : normal(n), uv(u), matIdx(idx) {}
	float3 normal = float3(0.0f); // 12 bytes
	float2 uv = float2(0.0f); // 8 bytes
	int matIdx; // 4 bytes
	// 24 bytes in total
};