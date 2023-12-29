#pragma once

#include "helper.h"
#include "texture.h"
#include "material.h"
#include "hit_info.h"
#include "base_scene.h"
#include "primitive_scene.h"
#include "file_scene.h"
#include "tlas_file_scene.h"

#define EPSILON	0.001f

namespace Tmpl8
{
	class Renderer : public TheApp
	{
	private:
		float3 m_camPositionToSet = 0;
		float3 m_camTargetToSet = 0;
		float m_avg = 10;
		float m_alpha = 1;
		float m_fps = 0;
		float m_rps = 0;
		bool m_inspectTraversal = false;
		float3 GetEdgeDebugColor(float2 uv);
	public:
		// game flow methods
		void Init();
		void ClearAccumulator();
		float3 HandleMirror(const Ray& ray, uint& seed, const float3& I, const float3& N, const int depth);
		float3 HandleDielectric(const Ray& ray, uint& seed, const float3& I, const float3& N, const int depth);
		float3 Sample(Ray& ray, uint& seed, int depth = 0);
		void ProcessTile(int tx, int ty, float& sum);
		void Tick( float deltaTime );
		void UI();
		void Shutdown() { /* implement if you want to do things on shutdown */ }
		// input handling
		void MouseUp( int button ) { /* implement if you want to detect mouse button presses */ }
		void MouseDown( int button ) { /* implement if you want to detect mouse button presses */ }
		void MouseMove( int x, int y ) { mousePos.x = x, mousePos.y = y; }
		void MouseWheel( float y ) { /* implement if you want to handle the mouse wheel */ }
		void KeyUp( int key ) { /* implement if you want to handle keys */ }
		void KeyDown( int key ) { /* implement if you want to handle keys */ }
		// data members
		int2 mousePos;
		float4* accumulator;
		TLASFileScene scene = TLASFileScene("../assets/scenes/inside_scene.xml");
		Camera camera;
		int spp = 1, passes = 1;
		bool animating = false;
		float energy, anim_time = 0;
		int depthLimit = 5;
	};
} // namespace Tmpl8