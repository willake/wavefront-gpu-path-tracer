#pragma once

#include "helper.h"
#include "texture.h"
#include "material.h"
#include "mesh.h"
#include "hit_info.h"
#include "base_scene.h"
#include "primitive_scene.h"
#include "tlas_file_scene.h"
#include "scene_buffer.h"

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
    uint m_extensionRayCount = 0;
    uint m_shadowRayCount = 0;
    bool m_disableUI = false;
    Buffer *GetPrimaryRayBuffer();
    Buffer *GetExtensionRayBuffer();
    void SwitchPrimaryRayBuffer();

  public:
    Renderer();
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
    Camera camera;
    int spp = 1, passes = 1;
    bool animating = false;
    float energy, anim_time = 0;
    int depthLimit = 5;

    Kernel *kernelGeneratePrimaryRays;
    Kernel *kernelExtend;
    Kernel *kernelShade;
    Kernel *kernelConnect;
    Kernel *kernelFinalize;
    Kernel *kernelClearAccumulator;
    Kernel *kernelDebugRenderTraversal; // debug

    float4 *Ts;
    float4 *Es;
    Ray *rays1;
    Ray *rays2; // have 2 ray buffers for switching, preventing writing the same buffer
    bool useRays1AsPrimary = true;
    ShadowRay *shadowrays;
    uint extensionCounter = 0;
    uint shadowrayCounter = 0;
    Buffer *rayBuffer1;
    Buffer *rayBuffer2;
    Buffer *extensionCounterBuffer;
    Buffer *shadowrayBuffer;
    Buffer *shadowrayCounterBuffer;
    Buffer *accumulatorBuffer;
    Buffer *seedBuffer;
    Buffer *screenPixelBuffer;
    Buffer *TBuffer;
    Buffer *EBuffer;
    Buffer *renderStateBuffer;

    SceneBuffer *sceneBuffer;
    string scenePath = "../assets/scenes/hollow_knight_scene.xml";
};
} // namespace Tmpl8