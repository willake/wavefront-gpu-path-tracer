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
    kernelConnect = new Kernel("../cl/kernels.cl", "connect");

    seeds = new uint[SCRWIDTH * SCRHEIGHT];
    pixels = new float4[SCRWIDTH * SCRHEIGHT];
    rays = new Ray[SCRWIDTH * SCRHEIGHT];
    rayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray), rays);
    // gpuaccumulator = new uint[SCRWIDTH * SCRHEIGHT];
    accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), accumulator);
    seedBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(uint), seeds);
    pixelBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), pixels);

    // init extension rays
    extensionrays = new Ray[SCRWIDTH * SCRHEIGHT];
    extensionrayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray), extensionrays);
    extensionCounterBuffer = new Buffer(sizeof(uint), &extensionCounter);

    // init shadow rays
    shadowrays = new ShadowRay[SCRWIDTH * SCRHEIGHT];
    shadowrayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(ShadowRay), shadowrayBuffer);
    shadowrayCounterBuffer = new Buffer(sizeof(uint), &shadowrayCounter);
}

void Renderer::ClearAccumulator()
{
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
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
    kernelGeneratePrimaryRays->SetArguments(rayBuffer, seedBuffer, pixelBuffer, SCRWIDTH, SCRHEIGHT, camera.camPos,
                                            camera.topLeft, camera.topRight, camera.bottomLeft, spp);
    rayBuffer->CopyToDevice(true);
    seedBuffer->CopyToDevice(true);
    pixelBuffer->CopyToDevice(true);
    kernelGeneratePrimaryRays->Run(SCRWIDTH * SCRHEIGHT);
    kernelExtend->SetArguments(rayBuffer, scene.triBuffer, scene.triIdxBuffer, scene.bvhNodeBuffer, scene.bvhBuffer,
                               scene.blasBuffer, scene.tlasNodeBuffer, scene.meshInsBuffer, scene.lightBuffer,
                               (int)scene.lightCount);
    kernelExtend->Run(SCRWIDTH * SCRHEIGHT);
    //// accumulatorBuffer->CopyToDevice(true);
    extensionCounter = 0;
    extensionCounterBuffer->CopyToDevice();
    shadowrayCounter = 0;
    shadowrayCounterBuffer->CopyToDevice();
    kernelShade->SetArguments(pixelBuffer, rayBuffer, seedBuffer, scene.skydomeBuffer, scene.skydome.width,
                              scene.skydome.height, scene.floorBuffer, scene.triExBuffer, scene.blasBuffer,
                              scene.materialBuffer, scene.texturePixelBuffer, scene.textureBuffer, scene.lightBuffer,
                              (int)scene.lightCount, extensionrayBuffer, shadowrayBuffer, extensionCounterBuffer,
                              shadowrayCounterBuffer);
    kernelShade->Run(SCRWIDTH * SCRHEIGHT);

    int pa = 0;
    // run extension rays and shadow rays
    while (pa < 10 || extensionCounter > 0 || shadowrayCounter > 0)
    {
        uint extensionCount = extensionCounter;
        uint shadowCount = shadowrayCounter;
        extensionCounter = 0;
        extensionCounterBuffer->CopyToDevice();
        shadowrayCounter = 0;
        shadowrayCounterBuffer->CopyToDevice();

        if (extensionCount > 0)
        {
            kernelExtend->SetArguments(extensionrayBuffer, scene.triBuffer, scene.triIdxBuffer, scene.bvhNodeBuffer,
                                       scene.bvhBuffer, scene.blasBuffer, scene.tlasNodeBuffer, scene.meshInsBuffer,
                                       scene.lightBuffer, (int)scene.lightCount);
            kernelExtend->Run(extensionCount);
        }

        if (shadowCount > 0)
        {
            kernelShade->SetArguments(pixelBuffer, rayBuffer, seedBuffer, scene.skydomeBuffer, scene.skydome.width,
                                      scene.skydome.height, scene.floorBuffer, scene.triExBuffer, scene.blasBuffer,
                                      scene.materialBuffer, scene.texturePixelBuffer, scene.textureBuffer,
                                      scene.lightBuffer, extensionrayBuffer, shadowrayBuffer, extensionCounterBuffer,
                                      shadowrayCounterBuffer);
            kernelShade->Run(shadowCount);
        }
        extensionCounterBuffer->CopyFromDevice();
        shadowrayCounterBuffer->CopyFromDevice();
        pa++;
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
    ImGui::Text("Energy: %fk", energy / 1000);
    ImGui::Text("RPS: %.1f Mrays/s", m_rps);
    ImGui::Text("Camera Pos: (%.2f, %.2f, %.2f)", camera.camPos.x, camera.camPos.y, camera.camPos.z);
    ImGui::Text("Camera Target: (%.2f, %.2f, %.2f)", camera.camTarget.x, camera.camTarget.y, camera.camTarget.z);
    // reset accumulator if changes have been made
    if (changed) ClearAccumulator();
}