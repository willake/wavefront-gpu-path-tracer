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
// - Can be evaluated at arbitrary time - for motion blur
// - Has some high-frequency details - for filtering
// Some speed tricks that severely affect maintainability
// are enclosed in #ifdef SPEEDTRIX / #endif. Mind these
// if you plan to alter the scene in any way.
// -----------------------------------------------------------

#define PLANE_X(o,i) { float t =-(ray.O.x+o) * ray.rD.x; if(t<ray.t&&t>0)ray.t=t,ray.objIdx=i;}
#define PLANE_Y(o,i) { float t =-(ray.O.y+o) * ray.rD.y; if(t<ray.t&&t>0)ray.t=t,ray.objIdx=i;}
#define PLANE_Z(o,i) { float t =-(ray.O.z+o) * ray.rD.z; if(t<ray.t&&t>0)ray.t=t,ray.objIdx=i;}

namespace Tmpl8 {

	// -----------------------------------------------------------
	// Scene class
	// We intersect this. The query is internally forwarded to the
	// list of primitives, so that the nearest hit can be returned.
	// For this hit (distance, obj id), we can query the normal and
	// albedo.
	// -----------------------------------------------------------
	class Scene
	{
	public:
		Scene()
		{
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
			SetTime(0);
			// Note: once we have triangle support we should get rid of the class
			// hierarchy: virtuals reduce performance somewhat.
		}
		void SetTime(float t)
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
		float3 GetLightPos() const
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
		float3 RandomPointOnLight(const float r0, const float r1) const
		{
#ifndef FOURLIGHTS
			// get a random position on the swinging quad
			const float size = quad.size;
			float3 corner1 = TransformPosition(float3(-size, 0, -size), quad.T);
			float3 corner2 = TransformPosition(float3(size, 0, -size), quad.T);
			float3 corner3 = TransformPosition(float3(-size, 0, size), quad.T);
			return corner1 + r0 * (corner2 - corner1) + r1 * (corner3 - corner1);
#else
			// select a random light and use that
			uint lightIdx = (uint)(r0 * 4);
			const Quad& q = quad[lightIdx];
			// renormalize r0 for reuse
			float stratum = lightIdx * 0.25f;
			float r2 = (r0 - stratum) / (1 - stratum);
			// get a random position on the selected quad
			const float size = q.size;
			float3 corner1 = TransformPosition(float3(-size, 0, -size), q.T);
			float3 corner2 = TransformPosition(float3(size, 0, -size), q.T);
			float3 corner3 = TransformPosition(float3(-size, 0, size), q.T);
			return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
#endif
		}
		float3 RandomPointOnLight(uint& seed) const
		{
			return RandomPointOnLight(RandomFloat(seed), RandomFloat(seed));
		}
		void GetLightQuad(float3& v0, float3& v1, float3& v2, float3& v3, const uint idx = 0)
		{
#ifndef FOURLIGHTS
			// return the four corners of the swinging quad, clockwise, for solid angle sampling
			const Quad& q = quad;
#else
			// return four corners of the specified light
			const Quad& q = quad[idx];
#endif
			const float size = q.size;
			v0 = TransformPosition(float3(-size, 0, size), q.T);
			v1 = TransformPosition(float3(size, 0, size), q.T);
			v2 = TransformPosition(float3(size, 0, -size), q.T);
			v3 = TransformPosition(float3(-size, 0, -size), q.T);
		}
		float3 GetLightColor() const
		{
			return float3(24, 24, 22);
		}
		float3 GetAreaLightColor() const
		{
#ifdef FOURLIGHTS
			return quad[0].GetAlbedo(float3(0)); // they're all the same color
#else
			return quad.GetAlbedo(float3(0));
#endif
		}
		float GetLightArea() const
		{
#ifdef FOURLIGHTS
			return sqrf(quad[0].size * 2); // all the same size
#else
			return sqrf(quad.size * 2);
#endif
		}
		constexpr float GetLightCount() const
		{
#ifdef FOURLIGHTS
			return 4; // what did you expect
#else
			return 1;
#endif
		}
		void FindNearest(Ray& ray) const
		{
			// room walls - ugly shortcut for more speed
			if (ray.D.x < 0) PLANE_X(3, 4) else PLANE_X(-2.99f, 5);
			if (ray.D.y < 0) PLANE_Y(1, 6) else PLANE_Y(-2, 7);
			if (ray.D.z < 0) PLANE_Z(3, 8) else PLANE_Z(-3.99f, 9);
			quad.Intersect(ray);
			sphere.Intersect(ray);
			sphere2.Intersect(ray);
			cube.Intersect(ray);
			torus.Intersect(ray);
		}
		bool IsOccluded(const Ray& ray) const
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
		float3 GetNormal(const int objIdx, const float3 I, const float3 wo) const
		{
			// we get the normal after finding the nearest intersection:
			// this way we prevent calculating it multiple times.
			if (objIdx == -1) return float3(0); // or perhaps we should just crash
			float3 N;
#ifdef FOURLIGHTS
			if (objIdx == 0) N = quad[0].GetNormal(I); // they're all oriented the same
#else
			if (objIdx == 0) N = quad.GetNormal(I);
#endif
			else if (objIdx == 1) N = sphere.GetNormal(I);
			else if (objIdx == 2) N = sphere2.GetNormal(I);
			else if (objIdx == 3) N = cube.GetNormal(I);
			else if (objIdx == 10) N = torus.GetNormal(I);
			else
			{
				// faster to handle the 6 planes without a call to GetNormal
				N = float3(0);
				N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
			}
			if (dot(N, wo) > 0) N = -N; // hit backside / inside
			return N;
		}
		float3 GetAlbedo(int objIdx, float3 I) const
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
		float GetReflectivity(int objIdx, float3 I) const
		{
			if (objIdx == 1 /* ball */) return 1;
			if (objIdx == 6 /* floor */) return 0.3f;
			return 0;
		}
		float GetRefractivity(int objIdx, float3 I) const
		{
			return (objIdx == 3 || objIdx == 10) ? 1.0f : 0.0f;
		}
		float3 GetAbsorption(int objIdx)
		{
			return objIdx == 3 ? float3(0.5f, 0, 0.5f) : float3(0);
		}
		__declspec(align(64)) // start a new cacheline here
			float animTime = 0;
#ifdef FOURLIGHTS
		Quad quad[4];
#else
		Quad quad, dummyQuad1, dummyQuad2, dummyQuad3;
#endif
		Sphere sphere;
		Sphere sphere2;
		Cube cube;
		Plane plane[6];
		Torus torus;
	};
}