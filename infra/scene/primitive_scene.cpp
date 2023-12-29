#include "precomp.h"
#include "primitive_scene.h"

PrimitiveScene::PrimitiveScene()
{
	errorMaterial.albedo = float3(255, 192, 203) / 255.f;
	// we store all primitives in one continuous buffer
#ifdef FOURLIGHTS
	for (int i = 0; i < 4; i++) quad[i] = Quad(0, 0.5f);	// 0: four light sources
#else
	quad = Quad(0, 1);									// 0: light source
#endif
	sphere = Sphere(1, float3(0), 0.6f);				// 1: bouncing ball
	sphere2 = Sphere(2, float3(0, 2.5f, -3.07f), 8);	// 2: rounded corners
	cube = Cube(3, float3(0), float3(1.15f));			// 3: cube
	plane[0] = Plane(4, float3(1, 0, 0), 3);			// 4: left wall
	plane[1] = Plane(5, float3(-1, 0, 0), 2.99f);		// 5: right wall
	plane[2] = Plane(6, float3(0, 1, 0), 1);			// 6: floor
	plane[3] = Plane(7, float3(0, -1, 0), 2);			// 7: ceiling
	plane[4] = Plane(8, float3(0, 0, 1), 3);			// 8: front wall
	plane[5] = Plane(9, float3(0, 0, -1), 3.99f);		// 9: back wall
	torus = Torus(10, 0.8f, 0.25f);						// 10: torus
	torus.T = mat4::Translate(-0.25f, 0, 2) * mat4::RotateX(PI / 4);
	torus.invT = torus.T.Inverted();
	materials[0].isLight = true; // 0: light source
	materials[1].reflectivity = 1.0f; // 1: bouncing ball
	materials[2] = Material(); // 2: rounded corners
	materials[3].refractivity = 1.0f; // 3: cube
	materials[3].absorption = float3(0.5f, 0, 0.5f);
	materials[4] = Material(true); // 4: left wall
	materials[5] = Material(true); // 5: right wall
	materials[6] = Material(true); // 6: floor
	materials[6].reflectivity = 0.3f;
	materials[7] = Material(); // 7: ceiling
	materials[8] = Material(); // 8: front wall
	materials[9] = Material(); // 9: back wall
	materials[10].refractivity = 1.0f; // 10: torus
	SetTime(0);
	// Note: once we have triangle support we should get rid of the class
	// hierarchy: virtuals reduce performance somewhat.
}

void PrimitiveScene::SetTime(float t)
{
	// default time for the scene is simply 0. Updating/ the time per frame 
	// enables animation. Updating it per ray can be used for motion blur.
	animTime = t;
#ifdef FOURLIGHTS
	// four light sources are stationary
	quad[0].T = mat4::Translate(-1, 1.5f, -1), quad[0].invT = quad[0].T.FastInvertedTransformNoScale();
	quad[1].T = mat4::Translate(1, 1.5f, -1), quad[1].invT = quad[1].T.FastInvertedTransformNoScale();
	quad[2].T = mat4::Translate(1, 1.5f, 1), quad[2].invT = quad[2].T.FastInvertedTransformNoScale();
	quad[3].T = mat4::Translate(-1, 1.5f, 1), quad[3].invT = quad[3].T.FastInvertedTransformNoScale();
#else
	// light source animation: swing
	mat4 M1base = mat4::Translate(float3(0, 2.6f, 2));
	mat4 M1 = M1base * mat4::RotateZ(sinf(animTime * 0.6f) * 0.1f) * mat4::Translate(float3(0, -0.9f, 0));
	quad.T = M1, quad.invT = M1.FastInvertedTransformNoScale();
#endif
	// cube animation: spin
	mat4 M2base = mat4::RotateX(PI / 4) * mat4::RotateZ(PI / 4);
	mat4 M2 = mat4::Translate(float3(1.8f, 0, 2.5f)) * mat4::RotateY(animTime * 0.5f) * M2base;
	cube.M = M2, cube.invM = M2.FastInvertedTransformNoScale();
	// sphere animation: bounce
	float tm = 1 - sqrf(fmodf(animTime, 2.0f) - 1);
	sphere.pos = float3(-1.8f, -0.4f + tm, 1);
}

float3 PrimitiveScene::GetLightPos() const
{
#ifndef FOURLIGHTS
	// light point position is the middle of the swinging quad
	float3 corner1 = TransformPosition(float3(-0.5f, 0, -0.5f), quad.T);
	float3 corner2 = TransformPosition(float3(0.5f, 0, 0.5f), quad.T);
	return (corner1 + corner2) * 0.5f - float3(0, 0.01f, 0);
#else
	// function is not valid when using four lights; we'll return the origin
	return float3(0);
#endif
}

float3 PrimitiveScene::GetSkyColor(const Ray& ray) const
{
	return float3(0);
}

float3 PrimitiveScene::GetLightColor() const
{
	return float3(24, 24, 22);
}

void PrimitiveScene::FindNearest(Ray& ray)
{
	// room walls - ugly shortcut for more speed
#ifdef SPEEDTRIX
	// prefetching
	const float3 spos = sphere.pos;
	const float3 ro = ray.O, rd = ray.D;
	// TODO: the room is actually just an AABB; use slab test
	static const __m128 x4min = _mm_setr_ps(3, 1, 3, 1e30f);
	static const __m128 x4max = _mm_setr_ps(-2.99f, -2, -3.99f, 1e30f);
	static const __m128 idmin = _mm_castsi128_ps(_mm_setr_epi32(4, 6, 8, -1));
	static const __m128 idmax = _mm_castsi128_ps(_mm_setr_epi32(5, 7, 9, -1));
	static const __m128 zero4 = _mm_setzero_ps();
	const __m128 selmask = _mm_cmpge_ps(ray.D4, zero4);
	const __m128i idx4 = _mm_castps_si128(_mm_blendv_ps(idmin, idmax, selmask));
	const __m128 x4 = _mm_blendv_ps(x4min, x4max, selmask);
	const __m128 d4 = _mm_sub_ps(zero4, _mm_mul_ps(_mm_add_ps(ray.O4, x4), ray.rD4));
	const __m128 mask4 = _mm_cmple_ps(d4, zero4);
	const __m128 t4 = _mm_blendv_ps(d4, _mm_set1_ps(1e34f), mask4);
	/* first: unconditional */  ray.t = t4.m128_f32[0], ray.objIdx = idx4.m128i_i32[0];
	if (t4.m128_f32[1] < ray.t) ray.t = t4.m128_f32[1], ray.objIdx = idx4.m128i_i32[1];
	if (t4.m128_f32[2] < ray.t) ray.t = t4.m128_f32[2], ray.objIdx = idx4.m128i_i32[2];
#else
	if (ray.D.x < 0) PLANE_X(3, 4) else PLANE_X(-2.99f, 5);
	if (ray.D.y < 0) PLANE_Y(1, 6) else PLANE_Y(-2, 7);
	if (ray.D.z < 0) PLANE_Z(3, 8) else PLANE_Z(-3.99f, 9);
#endif
#ifdef FOURLIGHTS
#ifdef SPEEDTRIX
	// efficient four-quad intersection by Jesse Vrooman
	const __m128 t = _mm_div_ps(_mm_add_ps(_mm_set1_ps(ray.O.y),
		_mm_set1_ps(-1.5)), _mm_xor_ps(_mm_set1_ps(ray.D.y), _mm_set1_ps(-0.0)));
	const __m128 Ix = _mm_add_ps(_mm_add_ps(_mm_set1_ps(ray.O.x),
		_mm_set_ps(1, -1, -1, 1)), _mm_mul_ps(t, _mm_set1_ps(ray.D.x)));
	const __m128 Iz = _mm_add_ps(_mm_add_ps(_mm_set1_ps(ray.O.z),
		_mm_set_ps(1, 1, -1, -1)), _mm_mul_ps(t, _mm_set1_ps(ray.D.z)));
	const static __m128 size = _mm_set1_ps(0.25f);
	const static __m128 nsize = _mm_xor_ps(_mm_set1_ps(0.25f), _mm_set1_ps(-0.0));
	const __m128 maskedT = _mm_and_ps(t, _mm_and_ps(
		_mm_and_ps(_mm_cmpgt_ps(Ix, nsize), _mm_cmplt_ps(Ix, size)),
		_mm_and_ps(_mm_cmpgt_ps(Iz, nsize), _mm_cmplt_ps(Iz, size))));
	if (maskedT.m128_f32[3] > 0) ray.t = maskedT.m128_f32[3], ray.objIdx = 0;
	if (maskedT.m128_f32[2] > 0) ray.t = maskedT.m128_f32[2], ray.objIdx = 0;
	if (maskedT.m128_f32[1] > 0) ray.t = maskedT.m128_f32[1], ray.objIdx = 0;
	if (maskedT.m128_f32[0] > 0) ray.t = maskedT.m128_f32[0], ray.objIdx = 0;
#else
	for (int i = 0; i < 4; i++) quad[i].Intersect(ray);
#endif
#else
	quad.Intersect(ray);
#endif
#ifdef SPEEDTRIX // hardcoded spheres, a bit faster this way but very ugly
	{
		// SIMD sphere intersection code by Jesse Vrooman
		const __m128 oc = _mm_sub_ps(ray.O4, sphere.pos4);
		const float b = _mm_dp_ps(oc, ray.D4, 0x71).m128_f32[0];
		const float d = b * b - (_mm_dp_ps(oc, oc, 0x71).m128_f32[0] - 0.36f);
		if (d > 0)
		{
			const float t = -b - sqrtf(d);
			const bool hit = t < ray.t && t > 0;
			if (hit) { ray.t = t, ray.objIdx = 1; }
		};
	}
	{
		// SIMD sphere intersection code by Jesse Vrooman
		const static __m128 s4 = _mm_setr_ps(0, 2.5f, -3.07f, 1);
		const __m128 oc = _mm_sub_ps(ray.O4, s4);
		const float b = _mm_dp_ps(oc, ray.D4, 0x71).m128_f32[0];
		const float d = b * b - (_mm_dp_ps(oc, oc, 0x71).m128_f32[0] - 64.0f);
		if (d > 0)
		{
			const float t = sqrtf(d) - b;
			const bool hit = t < ray.t && t > 0;
			if (hit) { ray.t = t, ray.objIdx = 2; }
		};
	}
#else
	sphere.Intersect(ray);
	sphere2.Intersect(ray);
#endif
	cube.Intersect(ray);
	torus.Intersect(ray);
}

bool PrimitiveScene::IsOccluded(const Ray& ray)
{
	if (cube.IsOccluded(ray)) return true;
#ifdef SPEEDTRIX
	const float3 oc = ray.O - sphere.pos;
	const float b = dot(oc, ray.D), c = dot(oc, oc) - (0.6f * 0.6f);
	const float d = b * b - c;
	if (d > 0)
	{
		const float t = -b - sqrtf(d);
		const bool hit = t < ray.t && t > 0;
		if (hit) return true;
	}
#else
	if (sphere.IsOccluded(ray)) return true;
#endif
#ifdef FOURLIGHTS
	for (int i = 0; i < 4; i++) if (quad[i].IsOccluded(ray)) return true;
#else
	if (quad.IsOccluded(ray)) return true;
#endif
	if (torus.IsOccluded(ray)) return true;
	return false; // skip planes and rounded corners
}

HitInfo PrimitiveScene::GetHitInfo(const Ray& ray, const float3 I)
{
	HitInfo hitInfo = HitInfo(float3(0), float2(0), &errorMaterial);
	switch (ray.objIdx)
	{
	case 0:
		hitInfo.normal = quad.GetNormal(I);
		hitInfo.uv = float2(0);
		break;
	case 1:
		hitInfo.normal = sphere.GetNormal(I);
		hitInfo.uv = float2(0);
		break;
	case 2:
		hitInfo.normal = sphere2.GetNormal(I);
		hitInfo.uv = float2(0);
		break;
	case 3:
		hitInfo.normal = cube.GetNormal(I);
		hitInfo.uv = float2(0);
		break;
	case 10:
		hitInfo.normal = torus.GetNormal(I);
		hitInfo.uv = float2(0);
		break;
	default:
		hitInfo.normal = float3(0);
		hitInfo.normal[(ray.objIdx - 4) / 2] = 1 - 2 * (float)(ray.objIdx & 1);
		break;
	}
	if (dot(hitInfo.normal, ray.D) > 0) hitInfo.normal = -hitInfo.normal;
	hitInfo.material = &materials[ray.objIdx];
	return hitInfo;
}
//
//float3 PrimitiveScene::GetNormal(const int objIdx, const float3 I, const float3 wo) const
//{
//	// we get the normal after finding the nearest intersection:
//	// this way we prevent calculating it multiple times.
//	if (objIdx == -1) return float3(0); // or perhaps we should just crash
//	float3 N;
//#ifdef FOURLIGHTS
//	if (objIdx == 0) N = quad[0].GetNormal(I); // they're all oriented the same
//#else
//	if (objIdx == 0) N = quad.GetNormal(I);
//#endif
//	else if (objIdx == 1) N = sphere.GetNormal(I);
//	else if (objIdx == 2) N = sphere2.GetNormal(I);
//	else if (objIdx == 3) N = cube.GetNormal(I);
//	else if (objIdx == 10) N = torus.GetNormal(I);
//	else
//	{
//		// faster to handle the 6 planes without a call to GetNormal
//		N = float3(0);
//		N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
//	}
//	if (dot(N, wo) > 0) N = -N; // hit backside / inside
//	return N;
//}

float3 PrimitiveScene::GetAlbedo(int objIdx, float3 I) const
{
	if (objIdx == -1) return float3(0); // or perhaps we should just crash
#ifdef FOURLIGHTS
	if (objIdx == 0) return quad[0].GetAlbedo(I); // they're all the same
#else
	if (objIdx == 0) return quad.GetAlbedo(I);
#endif
	if (objIdx == 1) return sphere.GetAlbedo(I);
	if (objIdx == 2) return sphere2.GetAlbedo(I);
	if (objIdx == 3) return cube.GetAlbedo(I);
	if (objIdx == 10) return torus.GetAlbedo(I);
	return plane[objIdx - 4].GetAlbedo(I);
	// once we have triangle support, we should pass objIdx and the bary-
	// centric coordinates of the hit, instead of the intersection location.
}

int PrimitiveScene::GetTriangleCount() const
{
	return 0;
}