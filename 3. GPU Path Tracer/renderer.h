#pragma once

#include "helper.h"
#include "texture.h"
#include "material.h"
#include "mesh.h"
#include "hit_info.h"
#include "base_scene.h"
#include "primitive_scene.h"
#include "tlas_file_scene.h"

#define EPSILON 0.001f

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

  public:
    // game flow methods
    void Init();
    void ClearAccumulator();
    void Tick(float deltaTime);
    void UI();
    void Shutdown()
    { /* implement if you want to do things on shutdown */
    }
    // input handling
    void MouseUp(int button)
    { /* implement if you want to detect mouse button presses */
    }
    void MouseDown(int button)
    { /* implement if you want to detect mouse button presses */
    }
    void MouseMove(int x, int y)
    {
        mousePos.x = x, mousePos.y = y;
    }
    void MouseWheel(float y)
    { /* implement if you want to handle the mouse wheel */
    }
    void KeyUp(int key)
    { /* implement if you want to handle keys */
    }
    void KeyDown(int key)
    { /* implement if you want to handle keys */
    }
    // data members
    int2 mousePos;
    float4 *accumulator;
    uint *seeds;
    TLASFileScene scene = TLASFileScene("../assets/scenes/base_scene.xml");
    Camera camera;
    int spp = 1, passes = 1;
    bool animating = false;
    float energy, anim_time = 0;
    int depthLimit = 5;

    Kernel *kernelTestRayStructSize;
    Kernel *kernelGeneratePrimaryRays;
    Kernel *kernelExtend;
    Kernel *kernelShade;
    Kernel *kernelConnect;

    float4 *pixels;
    Ray *rays;
    Ray *extensionrays;
    ShadowRay *shadowrays;
    uint extensionCounter = 0;
    uint shadowrayCounter = 0;
    Buffer *rayBuffer;
    Buffer *extensionrayBuffer;
    Buffer *extensionCounterBuffer;
    Buffer *shadowrayBuffer;
    Buffer *shadowrayCounterBuffer;
    Buffer *accumulatorBuffer;
    Buffer *seedBuffer;
    Buffer *pixelBuffer; // final pixels for a tick
};
} // namespace Tmpl8