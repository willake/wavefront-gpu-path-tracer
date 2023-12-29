#pragma once

namespace Tmpl8 {
	//__declspec(align(128)) 
	__declspec(align(128)) class Ray
	{
	public:
		Ray() = default;
		Ray(const Ray& ray)
		{
			O = ray.O, D = ray.D, t = ray.t; rD = ray.rD; objIdx = ray.objIdx; barycentric = ray.barycentric; triIdx = ray.triIdx;
			inside = ray.inside;
		}
		Ray(const float3 origin, const float3 direction, const float distance = 1e34f, const int idx = -1)
		{
			O = origin, D = direction, t = distance;
			// calculate reciprocal ray direction for triangles and AABBs
			rD = float3(1 / D.x, 1 / D.y, 1 / D.z);
			objIdx = idx;
		}
		float3 IntersectionPoint() const { return O + t * D; }
		// ray data
		union { struct { float3 O; float d0; }; __m128 O4; };
		union { struct { float3 D; float d1; }; __m128 D4; };
		union { struct { float3 rD; float d2; }; __m128 rD4; };
		float t = 1e34f; // 4 bytes
		float2 barycentric = float2(0); // 8 bytes
		int objIdx = -1; // 4 bytes
		int triIdx = -1; // 4 bytes
		int traversed = 0;
		int tested = 0;
		bool inside = false; // 1 bytes // true when in medium
	};
}