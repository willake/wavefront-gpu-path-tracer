#include "precomp.h"
#include "helper.h"
#include "renderer.h"

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
    // create fp32 rgb pixel buffer to render to
    accumulator = (float4 *)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
    ClearAccumulator();

    kernelTestRayStructSize = new Kernel("../cl/kernels.cl", "testRayStructSize");
    kernelGeneratePrimaryRays = new Kernel("../cl/kernels.cl", "generatePrimaryRays");
    kernelExtend = new Kernel("../cl/extend.cl", "extend");
    kernelShade = new Kernel("../cl/shade.cl", "shade");
    kernelConnect = new Kernel("../cl/connect.cl", "connect");

    seeds = new uint[SCRWIDTH * SCRHEIGHT];
    pixels = new float4[SCRWIDTH * SCRHEIGHT];
    // gpuaccumulator = new uint[SCRWIDTH * SCRHEIGHT];
    accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), accumulator);
    seedBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(uint), seeds);
    pixelBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), pixels);

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
}

void Renderer::ClearAccumulator()
{
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
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

void Renderer::SwitchPrimaryRay()
{
    useRays1AsPrimary = !useRays1AsPrimary;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
    // animation
    if (animating) scene.SetTime(anim_time += deltaTime * 0.002f), ClearAccumulator();
    // pixel loop
    Timer t;

    // GPGPU
    // primary ray execution
    Buffer *primaryRayBuffer = GetPrimaryRayBuffer();
    Buffer *extensionRayBuffer = GetExtensionRayBuffer();
    kernelGeneratePrimaryRays->SetArguments(primaryRayBuffer, seedBuffer, pixelBuffer, SCRWIDTH, SCRHEIGHT,
                                            camera.camPos, camera.topLeft, camera.topRight, camera.bottomLeft, spp);
    rayBuffer1->CopyToDevice(true);
    seedBuffer->CopyToDevice(true);
    pixelBuffer->CopyToDevice(true);
    kernelGeneratePrimaryRays->Run(SCRWIDTH * SCRHEIGHT);
    kernelExtend->SetArguments(primaryRayBuffer, scene.triBuffer, scene.triIdxBuffer, scene.bvhNodeBuffer,
                               scene.bvhBuffer, scene.blasBuffer, scene.tlasNodeBuffer, scene.meshInsBuffer,
                               scene.lightBuffer, (int)scene.lightCount);
    kernelExtend->Run(SCRWIDTH * SCRHEIGHT);
    //// accumulatorBuffer->CopyToDevice(true);
    extensionCounter = 0;
    extensionCounterBuffer->CopyToDevice(true);
    shadowrayCounter = 0;
    shadowrayCounterBuffer->CopyToDevice(true);
    kernelShade->SetArguments(pixelBuffer, primaryRayBuffer, seedBuffer, scene.skydomeBuffer, scene.skydome.width,
                              scene.skydome.height, scene.floorBuffer, scene.triExBuffer, scene.blasBuffer,
                              scene.materialBuffer, scene.texturePixelBuffer, scene.textureBuffer, scene.lightBuffer,
                              (int)scene.lightCount, extensionRayBuffer, shadowrayBuffer, extensionCounterBuffer,
                              shadowrayCounterBuffer, 0);
    kernelShade->Run(SCRWIDTH * SCRHEIGHT);

    extensionCounterBuffer->CopyFromDevice(true);
    shadowrayCounterBuffer->CopyFromDevice(true);
    m_extensionRayCount = extensionCounter;
    m_shadowRayCount = shadowrayCounter;

    SwitchPrimaryRay();

    int depth = 1;
    // run extension rays and shadow rays
    while (depth < 5 && (extensionCounter > 0 || shadowrayCounter > 0))
    {
        int extensionCount = extensionCounter;
        int shadowCount = shadowrayCounter;
        extensionCounter = 0;
        extensionCounterBuffer->CopyToDevice(true);
        shadowrayCounter = 0;
        shadowrayCounterBuffer->CopyToDevice(true);

        primaryRayBuffer = GetPrimaryRayBuffer();
        extensionRayBuffer = GetExtensionRayBuffer();

        if (extensionCount > 0)
        {
            kernelExtend->SetArguments(primaryRayBuffer, scene.triBuffer, scene.triIdxBuffer, scene.bvhNodeBuffer,
                                       scene.bvhBuffer, scene.blasBuffer, scene.tlasNodeBuffer, scene.meshInsBuffer,
                                       scene.lightBuffer, (int)scene.lightCount);
            kernelExtend->Run(extensionCount);
        }

        if (shadowCount > 0)
        {
            // kernelConnect->SetArguments(shadowrayBuffer, pixelBuffer, scene.triBuffer, scene.triIdxBuffer,
            //                             scene.bvhNodeBuffer, scene.bvhBuffer, scene.blasBuffer, scene.tlasNodeBuffer,
            //                             scene.meshInsBuffer, scene.lightBuffer, (int)scene.lightCount);
            // kernelConnect->Run(shadowCount);
        }

        kernelShade->SetArguments(pixelBuffer, primaryRayBuffer, seedBuffer, scene.skydomeBuffer, scene.skydome.width,
                                  scene.skydome.height, scene.floorBuffer, scene.triExBuffer, scene.blasBuffer,
                                  scene.materialBuffer, scene.texturePixelBuffer, scene.textureBuffer,
                                  scene.lightBuffer, extensionRayBuffer, shadowrayBuffer, extensionCounterBuffer,
                                  shadowrayCounterBuffer, depth);
        kernelShade->Run(extensionCount);
        extensionCounterBuffer->CopyFromDevice(true);
        shadowrayCounterBuffer->CopyFromDevice(true);
        depth++;
    }
    pixelBuffer->CopyFromDevice(true);
    //   accumulatorBuffer->CopyFromDevice(true);

    float scale = 1.0f / (spp + passes);
    for (int y = 0; y < SCRHEIGHT; y++)
        for (int x = 0; x < SCRWIDTH; x++)
        {
            accumulator[x + y * SCRWIDTH] += pixels[x + y * SCRWIDTH];
            float4 pixel = accumulator[x + y * SCRWIDTH] * scale;
            screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
        }

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
    ImGui::Checkbox("Inspect Traversal", &m_inspectTraversal);
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
    Ray r = camera.GetPrimaryRay((float)mousePos.x, (float)mousePos.y);
    scene.FindNearest(r);
    ImGui::Text("Object id: %i", r.objIdx);
    ImGui::Text("Triangle count: %i", scene.GetTriangleCount());
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