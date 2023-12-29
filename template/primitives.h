#pragma once

// -----------------------------------------------------------
// scene.h
// Simple test scene for ray tracing experiments. Goals:
// - Super-fast scene intersection
// - Easy interface: scene.FindNearest / IsOccluded
// - With normals and albedo: GetNormal / GetAlbedo
// - Area light source (animated), for light transport
// - Primitives can be hit from inside - for dielectrics
// - Can be extended with other primitives and/or a BVH
// - Optionally animated - for temporal experiments
// - Not everything is axis aligned - for cache experiments
// - Can be 
// d at arbitrary time - for motion blur
// - Has some high-frequency details - for filtering
// Some speed tricks that severely affect maintainability
// are enclosed in #ifdef SPEEDTRIX / #endif. Mind these
// if you plan to alter the scene in any way.
// -----------------------------------------------------------

#define SPEEDTRIX
//#define FOURLIGHTS

namespace Tmpl8 {
	// -----------------------------------------------------------
// Sphere primitive
// Basic sphere, with explicit support for rays that start
// inside it. Good candidate for a dielectric material.
// -----------------------------------------------------------
	class Sphere
	{
	public:
		Sphere() = default;
		Sphere(int idx, float3 p, float r) :
			pos(p), r2(r* r), invr(1 / r), objIdx(idx) {}
		void Intersect(Ray& ray) const
		{
#if 0
			const __m128 oc = _mm_sub_ps(ray.O4, _mm_set_ps(1, 1, this->pos.y, -1.8f));
			const float b = _mm_dp_ps(oc, ray.D4, 0xFF).m128_f32[0];
			const float c = _mm_dp_ps(oc, oc, 0xFF).m128_f32[0];
#else
			float3 oc = ray.O - this->pos;
			float b = dot(oc, ray.D);
			float c = dot(oc, oc) - this->r2;
#endif
			float t, d = b * b - c;
			if (d <= 0) return;
			d = sqrtf(d), t = -b - d;
			bool hit = t < ray.t && t > 0;
			if (hit)
			{
				ray.t = t, ray.objIdx = objIdx;
				return;
			}
			if (c > 0) return; // we're outside; safe to skip option 2
			t = d - b, hit = t < ray.t && t > 0;
			if (hit) ray.t = t, ray.objIdx = objIdx;
		}
		bool IsOccluded(const Ray& ray) const
		{
#if 0
			const __m128 oc = _mm_sub_ps(ray.O4, _mm_set_ps(1, 1, this->pos.y, -1.8f));
			const float b = _mm_dp_ps(oc, ray.D4, 0xFF).m128_f32[0];
			const float c = _mm_dp_ps(oc, oc, 0xFF).m128_f32[0];
#else
			float3 oc = ray.O - this->pos;
			float b = dot(oc, ray.D);
			float c = dot(oc, oc) - this->r2;
#endif
			float t, d = b * b - c;
			if (d <= 0) return false;
			d = sqrtf(d), t = -b - d;
			bool hit = t < ray.t && t > 0;
			return hit;
		}
		float3 GetNormal(const float3 I) const
		{
			return (I - this->pos) * invr;
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(0.93f);
		}
		union
		{
			float3 pos;
			__m128 pos4;
		};
		float r2 = 0, invr = 0;
		int objIdx = -1;
	};

	// -----------------------------------------------------------
	// Plane primitive
	// Basic infinite plane, defined by a normal and a distance
	// from the origin (in the direction of the normal).
	// -----------------------------------------------------------
	class Plane
	{
	public:
		Plane() = default;
		Plane(int idx, float3 normal, float dist, float textureOffset = 1.0f) : N(normal), d(dist), objIdx(idx) {
			invto = 1.f / textureOffset;
		}
		void Intersect(Ray& ray) const
		{
			float t = -(dot(ray.O, this->N) + this->d) / (dot(ray.D, this->N));
			if (t < ray.t && t > 0) ray.t = t, ray.objIdx = objIdx;
		}
		float3 GetNormal(const float3 I) const
		{
			return N;
		}
		float2 GetUV(const float3 I) const
		{
			if (N.y == 1)
			{
				float u = I.x;
				float v = I.z;

				u *= invto;
				v *= invto;

				u = u - floor(u);
				v = v - floor(v);

				// Return the UV coordinates
				return float2(u, v);
			}
			return float2(0);
		}
		float3 GetAlbedo(const float3 I) const
		{
			if (N.y == 1)
			{
				// floor albedo: checkerboard
				int ix = (int)(I.x * 2 + 96.01f);
				int iz = (int)(I.z * 2 + 96.01f);
				// add deliberate aliasing to two tile
				if (ix == 98 && iz == 98) ix = (int)(I.x * 32.01f), iz = (int)(I.z * 32.01f);
				if (ix == 94 && iz == 98) ix = (int)(I.x * 64.01f), iz = (int)(I.z * 64.01f);
				return float3(((ix + iz) & 1) ? 1 : 0.3f);
			}
			else if (N.z == -1)
			{
				// back wall: logo
				static Surface logo("../assets/logo.png");
				int ix = (int)((I.x + 4) * (128.0f / 8)), iy = (int)((2 - I.y) * (64.0f / 3));
				uint p = logo.pixels[(ix & 127) + (iy & 63) * 128];
				uint3 i3((p >> 16) & 255, (p >> 8) & 255, p & 255);
				return float3(i3) * (1.0f / 255.0f);
			}
			else if (N.x == 1)
			{
				// left wall: red
				static Surface red("../assets/red.png");
				int ix = (int)((I.z - 4) * (512.0f / 7)), iy = (int)((2 - I.y) * (512.0f / 3));
				uint p = red.pixels[(ix & 511) + (iy & 511) * 512];
				uint3 i3((p >> 16) & 255, (p >> 8) & 255, p & 255);
				return float3(i3) * (1.0f / 255.0f);
			}
			else if (N.x == -1)
			{
				// right wall: blue
				static Surface blue("../assets/blue.png");
				int ix = (int)((I.z - 4) * (512.0f / 7)), iy = (int)((2 - I.y) * (512.0f / 3));
				uint p = blue.pixels[(ix & 511) + (iy & 511) * 512];
				uint3 i3((p >> 16) & 255, (p >> 8) & 255, p & 255);
				return float3(i3) * (1.0f / 255.0f);
			}
			return float3(0.93f);
		}
		float3 N;
		float d;
		int objIdx = -1;
		float invto = 1.f; // inversed texture offeset
	};

	// -----------------------------------------------------------
	// Cube primitive
	// Oriented cube. Unsure if this will also work for rays that
	// start inside it; maybe not the best candidate for testing
	// dielectrics.
	// -----------------------------------------------------------
	class Cube
	{
	public:
		Cube() = default;
		Cube(int idx, float3 pos, float3 size, mat4 transform = mat4::Identity())
		{
			objIdx = idx;
			b[0] = float4(pos - 0.5f * size, 1);
			b[1] = float4(pos + 0.5f * size, 1);
			M = transform, invM = transform.FastInvertedTransformNoScale();
		}
		void Intersect(Ray& ray) const
		{
			// 'rotate' the cube by transforming the ray into object space
			// using the inverse of the cube transform.
#ifdef SPEEDTRIX
	// an AABB can be efficiently tested with SIMD - and matrix math is easy too. Profit!
			const __m128 a4 = ray.O4, b4 = ray.D4;
			__m128 v0 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[0])), v1 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[4]));
			__m128 v2 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[8])), v3 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[12]));
			_MM_TRANSPOSE4_PS(v0, v1, v2, v3);
			__m128 o4 = _mm_add_ps(_mm_add_ps(v0, v1), _mm_add_ps(v2, v3));
			__m128 v4 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[0])), v5 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[4]));
			__m128 v6 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[8])), v7 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[12]));
			_MM_TRANSPOSE4_PS(v4, v5, v6, v7);
			__m128 d4 = _mm_add_ps(_mm_add_ps(v4, v5), v6);
			__m128 rd4 = _mm_div_ps(_mm_set1_ps(1.0f), d4); // _mm_rcp_ps( d4 ); // reduced precision unsufficient?
			// AABB test
			__m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4, o4), rd4);
			__m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4, o4), rd4);
			__m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
			float tmax = min(vmax4.m128_f32[0], min(vmax4.m128_f32[1], vmax4.m128_f32[2]));
			float tmin = max(vmin4.m128_f32[0], max(vmin4.m128_f32[1], vmin4.m128_f32[2]));
			if (tmin < tmax) if (tmin > 0)
			{
				if (tmin < ray.t) ray.t = tmin, ray.objIdx = objIdx;
			}
			else if (tmax > 0)
			{
				if (tmax < ray.t) ray.t = tmax, ray.objIdx = objIdx;
			}
#else
			float3 O = TransformPosition(ray.O, invM);
			float3 D = TransformVector(ray.D, invM);
			float rDx = 1 / D.x, rDy = 1 / D.y, rDz = 1 / D.z;
			int signx = D.x < 0, signy = D.y < 0, signz = D.z < 0;
			float tmin = (b[signx].x - O.x) * rDx;
			float tmax = (b[1 - signx].x - O.x) * rDx;
			float tymin = (b[signy].y - O.y) * rDy;
			float tymax = (b[1 - signy].y - O.y) * rDy;
			if (tmin > tymax || tymin > tmax) return;
			tmin = max(tmin, tymin), tmax = min(tmax, tymax);
			float tzmin = (b[signz].z - O.z) * rDz;
			float tzmax = (b[1 - signz].z - O.z) * rDz;
			if (tmin > tzmax || tzmin > tmax) return;
			tmin = max(tmin, tzmin), tmax = min(tmax, tzmax);
			if (tmin > 0)
			{
				if (tmin < ray.t) ray.t = tmin, ray.objIdx = objIdx;
			}
			else if (tmax > 0)
			{
				if (tmax < ray.t) ray.t = tmax, ray.objIdx = objIdx;
			}
#endif
		}
		bool IsOccluded(const Ray& ray) const
		{
#ifdef SPEEDTRIX
			// an AABB can be efficiently tested with SIMD - and matrix math is easy too. Profit!
			const __m128 a4 = ray.O4, b4 = ray.D4;
			__m128 v0 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[0])), v1 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[4]));
			__m128 v2 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[8])), v3 = _mm_mul_ps(a4, _mm_load_ps(&invM.cell[12]));
			_MM_TRANSPOSE4_PS(v0, v1, v2, v3);
			__m128 o4 = _mm_add_ps(_mm_add_ps(v0, v1), _mm_add_ps(v2, v3));
			__m128 v4 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[0])), v5 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[4]));
			__m128 v6 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[8])), v7 = _mm_mul_ps(b4, _mm_load_ps(&invM.cell[12]));
			_MM_TRANSPOSE4_PS(v4, v5, v6, v7);
			__m128 d4 = _mm_add_ps(_mm_add_ps(v4, v5), v6);
			__m128 rd4 = _mm_div_ps(_mm_set1_ps(1.0f), d4); // _mm_rcp_ps( d4 ); // reduced precision unsufficient?
			// AABB test
			__m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4, o4), rd4);
			__m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4, o4), rd4);
			__m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
			float tmax = min(vmax4.m128_f32[0], min(vmax4.m128_f32[1], vmax4.m128_f32[2]));
			float tmin = max(vmin4.m128_f32[0], max(vmin4.m128_f32[1], vmin4.m128_f32[2]));
			return tmax > 0 && tmin < tmax && tmin < ray.t;
#else
			float3 O = TransformPosition_SSE(ray.O4, invM);
			float3 D = TransformVector_SSE(ray.D4, invM);
			float rDx = 1 / D.x, rDy = 1 / D.y, rDz = 1 / D.z;
			float t1 = (b[0].x - O.x) * rDx, t2 = (b[1].x - O.x) * rDx;
			float t3 = (b[0].y - O.y) * rDy, t4 = (b[1].y - O.y) * rDy;
			float t5 = (b[0].z - O.z) * rDz, t6 = (b[1].z - O.z) * rDz;
			float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
			float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
			return tmax > 0 && tmin < tmax && tmin < ray.t;
#endif
		}
		float3 GetNormal(const float3 I) const
		{
			// transform intersection point to object space
			float3 objI = TransformPosition(I, invM);
			// determine normal in object space
			float3 N = float3(-1, 0, 0);
			float d0 = fabs(objI.x - b[0].x), d1 = fabs(objI.x - b[1].x);
			float d2 = fabs(objI.y - b[0].y), d3 = fabs(objI.y - b[1].y);
			float d4 = fabs(objI.z - b[0].z), d5 = fabs(objI.z - b[1].z);
			float minDist = d0;
			if (d1 < minDist) minDist = d1, N.x = 1;
			if (d2 < minDist) minDist = d2, N = float3(0, -1, 0);
			if (d3 < minDist) minDist = d3, N = float3(0, 1, 0);
			if (d4 < minDist) minDist = d4, N = float3(0, 0, -1);
			if (d5 < minDist) minDist = d5, N = float3(0, 0, 1);
			// return normal in world space
			return TransformVector(N, M);
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(1, 1, 1);
		}
#ifdef SPEEDTRIX
		union { float4 b[2]; struct { __m128 bmin4, bmax4; }; };
#else
		float3 b[2];
#endif
		mat4 M, invM;
		int objIdx = -1;
	};

	// -----------------------------------------------------------
	// Quad primitive
	// Oriented quad, intended to be used as a light source.
	// -----------------------------------------------------------
	class Quad
	{
	public:
		Quad() = default;
		Quad(int idx, float s, mat4 transform = mat4::Identity())
		{
			objIdx = idx;
			size = s * 0.5f;
			T = transform, invT = transform.FastInvertedTransformNoScale();
		}
		void Intersect(Ray& ray) const
		{
			const float Oy = invT.cell[4] * ray.O.x + invT.cell[5] * ray.O.y + invT.cell[6] * ray.O.z + invT.cell[7];
			const float Dy = invT.cell[4] * ray.D.x + invT.cell[5] * ray.D.y + invT.cell[6] * ray.D.z;
			const float t = Oy / -Dy;
			if (t < ray.t && t > 0)
			{
				const float Ox = invT.cell[0] * ray.O.x + invT.cell[1] * ray.O.y + invT.cell[2] * ray.O.z + invT.cell[3];
				const float Oz = invT.cell[8] * ray.O.x + invT.cell[9] * ray.O.y + invT.cell[10] * ray.O.z + invT.cell[11];
				const float Dx = invT.cell[0] * ray.D.x + invT.cell[1] * ray.D.y + invT.cell[2] * ray.D.z;
				const float Dz = invT.cell[8] * ray.D.x + invT.cell[9] * ray.D.y + invT.cell[10] * ray.D.z;
				const float Ix = Ox + t * Dx, Iz = Oz + t * Dz;
				if (Ix > -size && Ix < size && Iz > -size && Iz < size)
					ray.t = t, ray.objIdx = objIdx;
			}
		}
		bool IsOccluded(const Ray& ray) const
		{
			const float Oy = invT.cell[4] * ray.O.x + invT.cell[5] * ray.O.y + invT.cell[6] * ray.O.z + invT.cell[7];
			const float Dy = invT.cell[4] * ray.D.x + invT.cell[5] * ray.D.y + invT.cell[6] * ray.D.z;
			const float t = Oy / -Dy;
			if (t < ray.t && t > 0)
			{
				const float Ox = invT.cell[0] * ray.O.x + invT.cell[1] * ray.O.y + invT.cell[2] * ray.O.z + invT.cell[3];
				const float Oz = invT.cell[8] * ray.O.x + invT.cell[9] * ray.O.y + invT.cell[10] * ray.O.z + invT.cell[11];
				const float Dx = invT.cell[0] * ray.D.x + invT.cell[1] * ray.D.y + invT.cell[2] * ray.D.z;
				const float Dz = invT.cell[8] * ray.D.x + invT.cell[9] * ray.D.y + invT.cell[10] * ray.D.z;
				const float Ix = Ox + t * Dx, Iz = Oz + t * Dz;
				return Ix > -size && Ix < size && Iz > -size && Iz < size;
			}
			return false;
		}
		float3 GetNormal(const float3 I) const
		{
			// TransformVector( float3( 0, -1, 0 ), T ) 
			return float3(-T.cell[1], -T.cell[5], -T.cell[9]);
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(10);
		}
		float size;
		mat4 T, invT;
		int objIdx = -1;
	};

	// -----------------------------------------------------------
	// Torus primitive - Inigo Quilez, ShaderToy 4sBGDy
	// -----------------------------------------------------------
	class Torus
	{
	public:
		Torus() = default;
		Torus::Torus(int idx, float a, float b) : objIdx(idx)
		{
			rc2 = a * a, rt2 = b * b;
			T = invT = mat4::Identity();
			r2 = sqrf(a + b);
		}
		void Intersect(Ray& ray) const
		{
			// via: https://www.shadertoy.com/view/4sBGDy
			float3 O = TransformPosition_SSE(ray.O4, invT);
			float3 D = TransformVector_SSE(ray.D4, invT);
			// extension rays need double precision for the quadratic solver!
			double po = 1, m = dot(O, O), k3 = dot(O, D), k32 = k3 * k3;
			// bounding sphere test
			double v = k32 - m + r2;
			if (v < 0) return;
			// setup torus intersection
			double k = (m - rt2 - rc2) * 0.5, k2 = k32 + rc2 * D.z * D.z + k;
			double k1 = k * k3 + rc2 * O.z * D.z, k0 = k * k + rc2 * O.z * O.z - rc2 * rt2;
			// solve quadratic equation
			if (fabs(k3 * (k32 - k2) + k1) < 0.0001)
			{
				swap(k1, k3);
				po = -1, k0 = 1 / k0, k1 = k1 * k0, k2 = k2 * k0, k3 = k3 * k0, k32 = k3 * k3;
			}
			double c2 = 2 * k2 - 3 * k32, c1 = k3 * (k32 - k2) + k1;
			double c0 = k3 * (k3 * (-3 * k32 + 4 * k2) - 8 * k1) + 4 * k0;
			c2 *= 0.33333333333, c1 *= 2, c0 *= 0.33333333333;
			double Q = c2 * c2 + c0, R = 3 * c0 * c2 - c2 * c2 * c2 - c1 * c1;
			double h = R * R - Q * Q * Q, z;
			if (h < 0)
			{
				const double sQ = sqrt(Q);
				z = 2 * sQ * cos(acos(R / (sQ * Q)) * 0.33333333333);
			}
			else
			{
#if 1
				const double sQ = cbrtFast(sqrt(h) + fabs(R)); // by Jeppe Vroegindeweij
#else
				const double sQ = cbrt(sqrt(h) + fabs(R)); // pow( sqrt( h ) + fabs( R ), 0.3333333 );
#endif
				z = copysign(fabs(sQ + Q / sQ), R);
			}
			z = c2 - z;
			double d1 = z - 3 * c2, d2 = z * z - 3 * c0;
			if (fabs(d1) < 1.0e-8)
			{
				if (d2 < 0) return;
				d2 = sqrt(d2);
			}
			else
			{
				if (d1 < 0) return;
				d1 = sqrt(d1 * 0.5), d2 = c1 / d1;
			}
			double t = 1e20;
			h = d1 * d1 - z + d2;
			if (h > 0)
			{
				h = sqrt(h);
				double t1 = -d1 - h - k3, t2 = -d1 + h - k3;
				t1 = (po < 0) ? 2 / t1 : t1, t2 = (po < 0) ? 2 / t2 : t2;
				if (t1 > 0) t = t1;
				if (t2 > 0) t = min(t, t2);
			}
			h = d1 * d1 - z - d2;
			if (h > 0)
			{
				h = sqrt(h);
				double t1 = d1 - h - k3, t2 = d1 + h - k3;
				t1 = (po < 0) ? 2 / t1 : t1, t2 = (po < 0) ? 2 / t2 : t2;
				if (t1 > 0) t = min(t, t1);
				if (t2 > 0) t = min(t, t2);
			}
			float ft = (float)t;
			if (ft > 0 && ft < ray.t) ray.t = ft, ray.objIdx = objIdx;
		}
		bool IsOccluded(const Ray& ray) const
		{
			// via: https://www.shadertoy.com/view/4sBGDy
			float3 O = TransformPosition_SSE(ray.O4, invT);
			float3 D = TransformVector_SSE(ray.D4, invT);
			float po = 1, m = dot(O, O), k3 = dot(O, D), k32 = k3 * k3;
			// bounding sphere test
			float v = k32 - m + r2;
			if (v < 0.0) return false;
			// setup torus intersection
			float k = (m - rt2 - rc2) * 0.5f, k2 = k32 + rc2 * D.z * D.z + k;
			float k1 = k * k3 + rc2 * O.z * D.z, k0 = k * k + rc2 * O.z * O.z - rc2 * rt2;
			// solve quadratic equation
			if (fabs(k3 * (k32 - k2) + k1) < 0.01f)
			{
				swap(k1, k3);
				po = -1, k0 = 1 / k0, k1 = k1 * k0, k2 = k2 * k0, k3 = k3 * k0, k32 = k3 * k3;
			}
			float c2 = 2 * k2 - 3 * k32, c1 = k3 * (k32 - k2) + k1;
			float c0 = k3 * (k3 * (-3 * k32 + 4 * k2) - 8 * k1) + 4 * k0;
			c2 *= 0.33333333333f, c1 *= 2, c0 *= 0.33333333333f;
			float Q = c2 * c2 + c0, R = 3 * c0 * c2 - c2 * c2 * c2 - c1 * c1;
			float h = R * R - Q * Q * Q, z = 0;
			if (h < 0)
			{
				const float sQ = sqrtf(Q);
				z = 2 * sQ * cosf(acosf(R / (sQ * Q)) * 0.3333333f);
			}
			else
			{
#if 1
				const float sQ = (float)cbrtFast(sqrtf(h) + fabsf(R));
#else
				const float sQ = cbrtf(sqrtf(h) + fabs(R)); // powf( sqrtf( h ) + fabs( R ), 0.3333333f );
#endif
				z = copysign(fabs(sQ + Q / sQ), R);
			}
			z = c2 - z;
			float d1 = z - 3 * c2, d2 = z * z - 3 * c0;
			if (fabs(d1) < 1.0e-4f)
			{
				if (d2 < 0) return false;
				d2 = sqrtf(d2);
			}
			else
			{
				if (d1 < 0.0) return false;
				d1 = sqrtf(d1 * 0.5f), d2 = c1 / d1;
			}
			float t = 1e20f;
			h = d1 * d1 - z + d2;
			if (h > 0)
			{
				float t1 = -d1 - sqrtf(h) - k3;
				t1 = (po < 0) ? 2 / t1 : t1;
				if (t1 > 0 && t1 < ray.t) return true;
			}
			h = d1 * d1 - z - d2;
			if (h > 0)
			{
				float t1 = d1 - sqrtf(h) - k3;
				t1 = (po < 0) ? 2 / t1 : t1;
				if (t1 > 0 && t1 < ray.t) return true;
			}
			return false;
		}
		float3 GetNormal(const float3 I) const
		{
			float3 L = TransformPosition(I, invT);
			float3 N = normalize(L * (dot(L, L) - rt2 - rc2 * float3(1, 1, -1)));
			return TransformVector(N, T);
		}
		float3 Torus::GetAlbedo(const float3 I) const
		{
			return float3(1); // material.albedo;
		}
		float rt2, rc2, r2;
		int objIdx;
		mat4 T, invT;
		// these are helper functions for the torus code
		// these function will find the cubic root up till a certain accuracy using the newtonian method
		float cbrtfFast(const float n) const {

			float x1 = n / 10.0f, x2 = 1.0f;
			int turn = 0;
			while (fabs(x1 - x2) > 0.00000001 && turn++ < 100)
				x1 = x2, x2 = (2.0f / 3.0f * x1) + (n / (3.0f * x1 * x1));
			return x2;
		}
		double cbrtFast(const double n) const {

			double x1 = n / 10.0f, x2 = 1.0f;
			int turn = 0;
			while (fabs(x1 - x2) > 0.00000001 && turn++ < 100)
				x1 = x2, x2 = (2.0f / 3.0f * x1) + (n / (3.0f * x1 * x1));
			return x2;
		}
	};
}