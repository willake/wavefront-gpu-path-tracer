#include "precomp.h"
#include "helper.h"
#include "renderer.h"

Renderer::Renderer()
{
    sceneBuffer = new SceneBuffer();
    TLASFileScene scene = TLASFileScene(scenePath, sceneBuffer);

    sceneBuffer->CreateBuffers();
    sceneBuffer->CopyToDevice();
}

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
    // create fp32 rgb pixel buffer to render to
    // accumulator = (float4 *)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
    accumulator = new float4[SCRWIDTH * SCRHEIGHT];

    kernelGeneratePrimaryRays = new Kernel("../cl/generate.cl", "generatePrimaryRays");
    kernelClearAccumulator = new Kernel("../cl/generate.cl", "clearAccumulator");
    kernelExtend = new Kernel("../cl/extend.cl", "extend");
    kernelShade = new Kernel("../cl/shade.cl", "shade");
    kernelConnect = new Kernel("../cl/connect.cl", "connect");
    kernelFinalize = new Kernel("../cl/finalize.cl", "finalize");
    kernelDebugRenderTraversal = new Kernel("../cl/debug.cl", "renderTraversal");

    seeds = new uint[SCRWIDTH * SCRHEIGHT];
    Ts = new float4[SCRWIDTH * SCRHEIGHT];
    Es = new float4[SCRWIDTH * SCRHEIGHT];
    // gpuaccumulator = new uint[SCRWIDTH * SCRHEIGHT];
    accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), accumulator);
    accumulatorBuffer->CopyToDevice(true);
    screenPixelBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(uint), screen->pixels);
    screenPixelBuffer->CopyToDevice(true);
    seedBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(uint), seeds);
    TBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), Ts);
    EBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), Es);
    TBuffer->CopyToDevice();
    EBuffer->CopyToDevice();

    // init extension rays
    rays1 = new Ray[SCRWIDTH * SCRHEIGHT];
    rayBuffer1 = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray), rays1);
    rays2 = new Ray[SCRWIDTH * SCRHEIGHT];
    rayBuffer2 = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray), rays2);
    extensionCounterBuffer = new Buffer(sizeof(uint), &extensionCounter);

    // init shadow rays
    shadowrays = new ShadowRay[SCRWIDTH * SCRHEIGHT];
    shadowrayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(ShadowRay), shadowrays);
    shadowrayCounterBuffer = new Buffer(sizeof(uint), &shadowrayCounter);

    ClearAccumulator();
}

void Renderer::ClearAccumulator()
{
    kernelClearAccumulator->SetArguments(accumulatorBuffer);
    kernelClearAccumulator->Run(SCRWIDTH * SCRHEIGHT);
    // memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
    // accumulatorBuffer->CopyToDevice(true);
}

Buffer *Renderer::GetPrimaryRayBuffer()
{
    if (useRays1AsPrimary) { return rayBuffer1; }
    else { return rayBuffer2; }
}

Buffer *Renderer::GetExtensionRayBuffer()
{
    if (useRays1AsPrimary) { return rayBuffer2; }
    else { return rayBuffer1; }
}

void Renderer::SwitchPrimaryRayBuffer()
{
    useRays1AsPrimary = !useRays1AsPrimary;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
    // animation
    // if (animating) sceneBuffer->SetTime(anim_time += deltaTime * 0.002f), ClearAccumulator();
    // pixel loop
    Timer t;

    // GPGPU
    // primary ray execution
    Buffer *primaryRayBuffer = GetPrimaryRayBuffer();
    Buffer *extensionRayBuffer = GetExtensionRayBuffer();
    kernelGeneratePrimaryRays->SetArguments(TBuffer, EBuffer, primaryRayBuffer, seedBuffer, SCRWIDTH, SCRHEIGHT,
                                            camera.camPos, camera.topLeft, camera.topRight, camera.bottomLeft, spp);
    rayBuffer1->CopyToDevice(true);
    seedBuffer->CopyToDevice(true);
    kernelGeneratePrimaryRays->Run(SCRWIDTH * SCRHEIGHT);
    kernelExtend->SetArguments(primaryRayBuffer, sceneBuffer->triBuffer, sceneBuffer->triIdxBuffer,
                               sceneBuffer->bvhNodeBuffer, sceneBuffer->bvhBuffer, sceneBuffer->blasBuffer,
                               sceneBuffer->tlasNodeBuffer, sceneBuffer->meshInsBuffer, sceneBuffer->lightBuffer,
                               (int)sceneBuffer->lightCount);
    kernelExtend->Run(SCRWIDTH * SCRHEIGHT);

    if (m_inspectTraversal)
    {
        kernelDebugRenderTraversal->SetArguments(EBuffer, primaryRayBuffer);
        kernelDebugRenderTraversal->Run(SCRWIDTH * SCRHEIGHT);
    }
    else
    {
        extensionCounter = SCRWIDTH * SCRHEIGHT;
        shadowrayCounter = 0;
        extensionCounterBuffer->CopyToDevice(true);
        shadowrayCounterBuffer->CopyToDevice(true);

        int depth = 0;
        // run extension rays and shadow rays
        while (extensionCounter > 0)
        {
            int extensionCount = extensionCounter;

            extensionCounter = 0;
            shadowrayCounter = 0;
            extensionCounterBuffer->CopyToDevice(true);
            shadowrayCounterBuffer->CopyToDevice(true);

            primaryRayBuffer = GetPrimaryRayBuffer();
            extensionRayBuffer = GetExtensionRayBuffer();

            kernelShade->SetArguments(TBuffer, EBuffer, primaryRayBuffer, seedBuffer, sceneBuffer->skydomeBuffer,
                                      (int)sceneBuffer->skydomeTexture.width, (int)sceneBuffer->skydomeTexture.height,
                                      sceneBuffer->floorBuffer, sceneBuffer->triExBuffer, sceneBuffer->blasBuffer,
                                      sceneBuffer->materialBuffer, sceneBuffer->texturePixelBuffer,
                                      sceneBuffer->textureBuffer, sceneBuffer->lightBuffer,
                                      (int)sceneBuffer->lightCount, extensionRayBuffer, shadowrayBuffer,
                                      extensionCounterBuffer, shadowrayCounterBuffer, (int)depth);
            kernelShade->Run(extensionCount);

            extensionCounterBuffer->CopyFromDevice(true);
            shadowrayCounterBuffer->CopyFromDevice(true);
            m_extensionRayCount = extensionCounter;
            m_shadowRayCount = shadowrayCounter;

            if (extensionCounter > 0)
            {
                kernelExtend->SetArguments(extensionRayBuffer, sceneBuffer->triBuffer, sceneBuffer->triIdxBuffer,
                                           sceneBuffer->bvhNodeBuffer, sceneBuffer->bvhBuffer, sceneBuffer->blasBuffer,
                                           sceneBuffer->tlasNodeBuffer, sceneBuffer->meshInsBuffer,
                                           sceneBuffer->lightBuffer, (int)sceneBuffer->lightCount);
                kernelExtend->Run(extensionCounter);
            }

            if (shadowrayCounter > 0)
            {
                kernelConnect->SetArguments(TBuffer, EBuffer, shadowrayBuffer, sceneBuffer->triBuffer,
                                            sceneBuffer->triIdxBuffer, sceneBuffer->bvhNodeBuffer,
                                            sceneBuffer->bvhBuffer, sceneBuffer->blasBuffer,
                                            sceneBuffer->tlasNodeBuffer, sceneBuffer->meshInsBuffer,
                                            sceneBuffer->lightBuffer, (int)sceneBuffer->lightCount);
                kernelConnect->Run(shadowrayCounter);
            }

            SwitchPrimaryRayBuffer();
            depth++;
        }
    }
    EBuffer->CopyFromDevice(true);
    //   accumulatorBuffer->CopyFromDevice(true);

    float scale = 1.0f / (spp + passes);
    kernelFinalize->SetArguments(screenPixelBuffer, accumulatorBuffer, EBuffer, scale, m_inspectTraversal);
    kernelFinalize->Run(SCRWIDTH * SCRHEIGHT);
    accumulatorBuffer->CopyFromDevice(true);
    screenPixelBuffer->CopyFromDevice(true);
    // for (int y = 0; y < SCRHEIGHT; y++)
    //     for (int x = 0; x < SCRWIDTH; x++)
    //     {
    //         accumulator[x + y * SCRWIDTH] += pixels[x + y * SCRWIDTH];
    //         float4 pixel = accumulator[x + y * SCRWIDTH] * scale;
    //         pixel = float4(pow(pixel.x, 1.0f / 2.2f), pow(pixel.y, 1.0f / 2.2f), pow(pixel.z, 1.0f / 2.2f),
    //                        0); // gamma correction
    //         screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
    //     }

    //// performance report - running average - ms, MRays/s
    m_avg = (1 - m_alpha) * m_avg + m_alpha * t.elapsed() * 1000;
    if (m_alpha > 0.05f) m_alpha *= 0.75f;
    m_fps = 1000.0f / m_avg, m_rps = (SCRWIDTH * SCRHEIGHT) / m_avg;
    //   handle user input
    if (camera.HandleInput(deltaTime)) { ClearAccumulator(); }
    else spp += passes;
}

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
    // animation toggle
    bool changed = ImGui::Checkbox("Animate scene", &animating);
    changed |= ImGui::Checkbox("Inspect Traversal", &m_inspectTraversal);
    changed |= ImGui::SliderInt("spp", &passes, 1, 4, "%i");
    ImGui::SliderFloat("Camera move speed", &camera.moveSpeed, 1.0f, 10.0f, "%.2f");
    ImGui::SliderFloat("Camera turn speed", &camera.turnSpeed, 1.0f, 10.0f, "%.2f");
    // camera position field
    ImGuiFloat3("Position", m_camPositionToSet);
    ImGuiFloat3("Rotation", m_camTargetToSet);
    if (ImGui::Button("Set Camera"))
    {
        // Button was clicked, perform action (e.g., reset values)
        camera.SetCameraState(m_camPositionToSet, m_camTargetToSet);
        ClearAccumulator();
    }
    // ray query on mouse
    // Ray r = camera.GetPrimaryRay((float)mousePos.x, (float)mousePos.y);
    // sceneBuffer->FindNearest(r);
    ImGui::Text("Object id: %i", 0);
    ImGui::Text("Triangle count: %i", sceneBuffer->totalTriangleCount);
    ImGui::Text("Frame: %5.2f ms (%.1ffps)", m_avg, m_fps);
    ImGui::Text("spp: %i", spp);
    // ImGui::Text("Energy: %fk", energy / 1000);
    ImGui::Text("Extension ray count: %i", m_extensionRayCount);
    ImGui::Text("Shadow ray count: %i", m_shadowRayCount);
    ImGui::Text("RPS: %.1f Mrays/s", m_rps);
    ImGui::Text("Camera Pos: (%.2f, %.2f, %.2f)", camera.camPos.x, camera.camPos.y, camera.camPos.z);
    ImGui::Text("Camera Target: (%.2f, %.2f, %.2f)", camera.camTarget.x, camera.camTarget.y, camera.camTarget.z);
    // reset accumulator if changes have been made
    if (changed) ClearAccumulator();
}