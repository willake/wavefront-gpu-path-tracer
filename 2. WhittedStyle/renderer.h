#pragma once

#include "helper.h"
#include "texture.h"
#include "material.h"
#include "mesh.h"
#include "hit_info.h"
#include "base_scene.h"
#include "primitive_scene.h"
#include "file_scene.h"
#include "tlas_file_scene.h"

#define EPSILON	0.001f
#define INVSCRWIDTH 1.0f / SCRWIDTH
#define INVSCRHEIGHT 1.0f / SCRHEIGHT

namespace Tmpl8
{
	struct Test
	{
		int sizeO;
		int sizeD;
		int sizerD;
		int sizet;
		int sizebary;
		int sizeObjIdx;
		int sizeTriIdx;
		int sizeInside;
		int sizePadding;
		int sizeTotal;
	};


	class Renderer : public TheApp
	{
	private:
		float3 m_camPositionToSet = 0;
		float3 m_camTargetToSet = 0;
		float m_avg = 10;
		float m_alpha = 1;
		float m_fps = 0;
		float m_rps = 0;
		float m_rayHitCount = 0;
		bool m_inspectTraversal = false;
		bool m_inspectIntersectionTest = false;
		float m_totalTraversal = 0;
		float m_totalTests = 0;
		float m_averageTraversal = 0;
		float m_averageTests = 0;
		float m_peakTraversal = 0;
		float m_peakTests = 0;
		// about the scene
		std::chrono::microseconds m_buildTime;
		uint m_triangleCount = 0;
		uint m_maxTreeDepth = 0;
		float3 GetEdgeDebugColor(float2 uv);
	public:
		// game flow methods
		void Init();
		float3 Trace(Ray& ray, int depth);
		float3 DirectIllumination(float3 I, float3 N);
		void Tick(float deltaTime);
		void UI();
		void Shutdown() { /* implement if you want to do things on shutdown */ }
		// input handling
		void MouseUp(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseDown(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
		void MouseWheel(float y) { /* implement if you want to handle the mouse wheel */ }
		void KeyUp(int key) { /* implement if you want to handle keys */ }
		void KeyDown(int key) { /* implement if you want to handle keys */ }
		// data members
		int2 mousePos;
		float4* accumulator;
		TLASFileScene scene = TLASFileScene("../assets/scenes/base_scene.xml");
		Camera camera;
		bool animating = false;
		float anim_time = 0;
		int depthLimit = 5;
		Kernel* kernelTestRayStructSize;
		Kernel* kernelGeneratePrimaryRays;
		Kernel* kernelExtend;
		Kernel* kernelShade;
		Kernel* kernelConnect;
		Kernel* kernelFinalize;
		Buffer* bvhBuffer;
		Buffer* rayBuffer;
		Test* tests;
		Ray* rays;
	};
} // namespace Tmpl8