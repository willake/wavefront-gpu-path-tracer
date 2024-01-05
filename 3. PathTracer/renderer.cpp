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

    rays = new Ray[SCRWIDTH * SCRHEIGHT];
    rayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray), rays);
    gpuaccumulator = new uint[SCRWIDTH * SCRHEIGHT];
    accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(uint), gpuaccumulator);
}

void Renderer::ClearAccumulator()
{
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
}

float3 Renderer::HandleMirror(const Ray &ray, uint &seed, const float3 &I, const float3 &N, const int depth)
{
    float3 R = reflect(ray.D, N);
    Ray r(I + R * EPSILON, R);
    return Sample(r, seed, depth + 1);
}

float3 Renderer::HandleDielectric(const Ray &ray, uint &seed, const float3 &I, const float3 &N, const int depth)
{
    float3 R = reflect(ray.D, N);
    Ray r(I + R * EPSILON, R);
    float n1 = ray.inside ? 1.2f : 1, n2 = ray.inside ? 1 : 1.2f;
    float eta = n1 / n2, cosi = dot(-ray.D, N);
    float cost2 = 1.0f - eta * eta * (1 - cosi * cosi);
    float Fr = 1;
    if (cost2 > 0)
    {
        float a = n1 - n2, b = n1 + n2, R0 = (a * a) / (b * b), c = 1 - cosi;
        Fr = R0 + (1 - R0) * (c * c * c * c * c);
        float3 T = eta * ray.D + ((eta * cosi - sqrtf(fabs(cost2))) * N);
        Ray t(I + T * EPSILON, T);
        t.inside = !ray.inside;
        if (RandomFloat(seed) > Fr)
            return Sample(t, seed, depth + 1);
    }
    return Sample(r, seed, depth + 1);
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Sample(Ray &ray, uint &seed, int depth)
{
    scene.FindNearest(ray);
    // if (ray.objIdx == -1) return float3(0);
    if (ray.objIdx == -1)
        return scene.GetSkyColor(ray); // or a fancy sky color
    if (depth >= depthLimit)
        return float3(0);
    float3 I = ray.O + ray.t * ray.D;
    HitInfo hitInfo = scene.GetHitInfo(ray, I);
    float3 N = hitInfo.normal;
    float2 uv = hitInfo.uv;
    Material *material = hitInfo.material;
    float3 albedo = material->isAlbedoOverridden ? scene.GetAlbedo(ray.objIdx, I) : material->GetAlbedo(uv);

    /* visualize edges */    // return GetEdgeDebugColor(ray.barycentric);
    /* visualize normal */  // return (N + 1) * 0.5f;
    /* visualize distance */ // return 0.1f * float3( ray.t, ray.t, ray.t );
    /* visualize albedo */   // return albedo;
    // if (m_inspectTraversal) return GetTraverseCountColor(ray.traversed, );

    if (material->isLight)
        return scene.GetLightColor();

    float3 out_radiance(0);
    float reflectivity = material->reflectivity;
    float refractivity = material->refractivity;
    // float diffuseness = 1 - (reflectivity + refractivity);

    float3 medium_scale(1);
    if (ray.inside)
    {
        float3 absorption = material->absorption;
        medium_scale = expf(absorption * -ray.t);
    }

    // choose a type of transport
    float r = RandomFloat(seed);
    if (r < reflectivity) // handle pure speculars
    {
        return albedo * medium_scale * HandleMirror(ray, seed, I, N, depth);
    }
    else if (r < reflectivity + refractivity) // handle dielectrics
    {
        return albedo * medium_scale * HandleDielectric(ray, seed, I, N, depth);
    }
    else // diffuse surface
    {
        float3 R = diffusereflection(N, seed);
        float3 brdf = albedo * INVPI;
        Ray r(I + R * EPSILON, R);
        return medium_scale * brdf * 2 * PI * dot(R, N) * Sample(r, seed, depth + 1);
    }
}

float3 Renderer::GetEdgeDebugColor(float2 uv)
{
    if (abs(uv.x) < 0.03f || abs(uv.x - 1) < 0.03f || abs(uv.y) < 0.03f || abs(uv.y - 1) < 0.03f)
    {
        return float3(0, 0, 0);
    }
    else
    {
        return float3(1);
    }
}

// -----------------------------------------------------------
// Draw an 16x16 tile of pixels
// -----------------------------------------------------------
void Renderer::ProcessTile(int tx, int ty, float &sum)
{
    float scale = 1.0f / (spp + passes);
    uint seed = InitSeed(tx + ty * SCRWIDTH + spp * 1799);
    for (int y = ty * 16, v = 0; v < 16; v++, y++)
        for (int x = tx * 16, u = 0; u < 16; u++, x++)
        {
            for (int p = 0; p < passes; p++)
                accumulator[x + y * SCRWIDTH] += float4(
                    Sample(camera.GetPrimaryRay((float)x + RandomFloat(seed), (float)y + RandomFloat(seed)), seed), 0);
            float4 pixel = accumulator[x + y * SCRWIDTH] * scale;
            sum += pixel.x + pixel.y + pixel.z;
            screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
        }
}
static struct TileJob : public Job
{
    void Main()
    {
        renderer->ProcessTile(tx, ty, sum);
    }
    Job *Init(Renderer *r, int x, int y)
    {
        renderer = r, tx = x, ty = y, sum = 0;
        return this;
    }
    Renderer *renderer;
    int tx, ty;
    float sum;
} tileJob[4096];

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
    // animation
    if (animating)
        scene.SetTime(anim_time += deltaTime * 0.002f), ClearAccumulator();
    // pixel loop
    Timer t;

    // GPGPU
    kernelGeneratePrimaryRays->SetArguments(rayBuffer, SCRWIDTH, SCRHEIGHT, camera.camPos, camera.topLeft,
                                            camera.topRight, camera.bottomLeft);
    rayBuffer->CopyToDevice(true);
    kernelGeneratePrimaryRays->Run(SCRWIDTH * SCRHEIGHT);
    kernelExtend->SetArguments(rayBuffer, scene.triBuffer, scene.triIdxBuffer, scene.bvhNodeBuffer, scene.bvhBuffer,
                               scene.blasBuffer, scene.tlasNodeBuffer);
    kernelExtend->Run(SCRWIDTH * SCRHEIGHT);
    accumulatorBuffer->CopyToDevice(true);
    kernelShade->SetArguments(accumulatorBuffer, rayBuffer, scene.skydomeBuffer, scene.skydome.width,
                              scene.skydome.height, scene.floorBuffer, scene.triExBuffer, scene.blasBuffer);
    kernelShade->Run(SCRWIDTH * SCRHEIGHT);
    accumulatorBuffer->CopyFromDevice(true);

    for (int y = 0; y < SCRHEIGHT; y++)
        for (int x = 0; x < SCRWIDTH; x++)
        {
            screen->pixels[x + y * SCRWIDTH] = gpuaccumulator[x + y * SCRWIDTH];
        }

    // kernelExtend->Run(SCRWIDTH * SCRHEIGHT);
    // kernelConnect->Run(SCRWIDTH * SCRHEIGHT);
    // rayBuffer->CopyFromDevice(true);

    // render tiles using a simple job system
    // for (int jobIdx = 0, y = 0; y < SCRHEIGHT / 16; y++) for (int x = 0; x < SCRWIDTH / 16; x++)
    //	jm->AddJob2(tileJob[jobIdx++].Init(this, x, y));
    // jm->RunJobs();
    //// gather energy received on tiles
    // energy = 0;
    // for (int tiles = (SCRWIDTH / 16) * (SCRHEIGHT / 16), i = 0; i < tiles; i++)
    //	energy += tileJob[i].sum;
    //// performance report - running average - ms, MRays/s
    // m_avg = (1 - m_alpha) * m_avg + m_alpha * t.elapsed() * 1000;
    // if (m_alpha > 0.05f) m_alpha *= 0.75f;
    // m_fps = 1000.0f / m_avg, m_rps = (SCRWIDTH * SCRHEIGHT) / m_avg;
    //  handle user input
    if (camera.HandleInput(deltaTime))
    {
        ClearAccumulator();
    }
    else
        spp += passes;
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
    if (changed)
        ClearAccumulator();
}