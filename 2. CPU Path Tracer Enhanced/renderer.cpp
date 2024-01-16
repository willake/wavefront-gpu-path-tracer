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
}

void Renderer::ClearAccumulator()
{
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
    spp = 1;
}

float3 Renderer::CalculateMicrofacetBRDF(float3 L, float3 V, float3 N, float metallic, float roughness,
                                         float3 baseColor, float reflectance)
{
    float3 H = normalize(V + L); // H = normalize(V + L);
    float NdotH = dot(N, H), NdotV = dot(N, V), NdotL = dot(N, L), VdotH = dot(V, H);

    float3 f0 = float3(0.16 * (reflectance * reflectance));
    f0 = mix(f0, baseColor, metallic);

    float3 F = FresnelSchlick(VdotH, f0);
    float D = DistributionGGX(NdotH, roughness);
    float G = GeometrySmith(NdotV, NdotL, roughness);

    float3 spec = F * G * D / 4.0f * NdotL * NdotV;

    baseColor *= float3(1.0) - F;

    baseColor *= (1.0 - metallic);

    float3 diffuse = baseColor * INVPI;

    return diffuse + spec;
}

void Renderer::NEE(uint &seed, const float3 &I, const float3 &V, const float3 &N, float3 &albedo, float3 &E, float3 &T)
{
    uint lightIdx;
    float3 randomLightPos = scene.RandomPointOnLight(seed, lightIdx);
    Light light = scene.GetLightByLightIdx(lightIdx);

    // Light related information
    float3 L = randomLightPos - I;
    float dist = length(L);
    L *= 1 / dist;
    float ndotl = dot(N, L);
    float nldotl = dot(light.GetNormal(I), -L);
    float A = light.size * light.size;
    Ray shadowRay = Ray(I + L * EPSILON, L, dist - 2 * EPSILON);
    if (ndotl > 0 && nldotl > 0)
    {
        if (!scene.IsOccluded(shadowRay))
        {
            float solidAngle = (nldotl * A) / (dist * dist);
            float3 mBRDF = CalculateMicrofacetBRDF(L, V, N, 0.1f, 1, albedo, 0);
            E += T * light.color * solidAngle * mBRDF * ndotl * scene.lightCount;
        }
    }
}

Ray Renderer::HandleMirror(const Ray &ray, uint &seed, const float3 &I, const float3 &N)
{
    float3 R = reflect(ray.D, N);
    Ray r(I + R * EPSILON, R);
    return r;
}

Ray Renderer::HandleDielectric(const Ray &ray, uint &seed, const float3 &I, const float3 &N)
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
        if (RandomFloat(seed) > Fr) return t;
    }
    return r;
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Sample(Ray &ray, uint &seed)
{
    float3 T = float3(1, 1, 1), E = float3(0, 0, 0);
    int depth = 0;
    bool lastSpecular = false;
    while (true)
    {
        scene.FindNearest(ray);
        if (ray.objIdx == -1) // hit skybox
        {
            E += T * scene.GetSkyColor(ray);
            break;
        }
        float3 I = ray.O + ray.t * ray.D;
        HitInfo hitInfo = scene.GetHitInfo(ray, I);
        float3 N = hitInfo.normal;
        float2 uv = hitInfo.uv;
        Material *material = hitInfo.material;
        float3 albedo = material->isAlbedoOverridden ? scene.GetAlbedo(ray.objIdx, I) : material->GetAlbedo(uv);
        float3 brdf = albedo * INVPI;

        // Part of NEE
        // return black if it is a light soucre
        if (material->isLight)
        {
            if (lastSpecular || (depth == 0 && dot(-ray.D, N) > 0))
            {
                E += T * scene.GetLightByObjIdx(ray.objIdx).color;
            }
            else E += T * float3(0);
            break;
        }

        float reflectivity = material->reflectivity;
        float refractivity = material->refractivity;

        float3 medium_scale(1);
        // beer's law
        if (ray.inside)
        {
            float3 absorption = material->absorption;
            medium_scale = expf(absorption * -ray.t);
        }

        // might modify E
        NEE(seed, I, -ray.D, N, albedo, E, T);

        float p = SurvivalProb(T);

        if (depth > 1 && p < RandomFloat(seed)) break;

        // choose a type of transport
        float r = RandomFloat(seed);
        if (r < reflectivity) // handle pure speculars
        {
            T *= albedo * medium_scale / p;
            ray = HandleMirror(ray, seed, I, N);
            lastSpecular = true;
        }
        else if (r < reflectivity + refractivity) // handle dielectrics
        {
            T *= albedo * medium_scale / p;
            ray = HandleDielectric(ray, seed, I, N);
            lastSpecular = true;
        }
        else // diffuse surface
        {
            float3 R = cosineweighteddiffusereflection(N, seed);
            float PDF = dot(N, R) / PI;
            float3 mBRDF = CalculateMicrofacetBRDF(R, -ray.D, N, 0.1f, 1, albedo, 0);
            ray = Ray(I + R * EPSILON, R);
            T *= medium_scale * mBRDF * dot(R, N) / PDF / p;
        }
        depth++;
    }

    return E;
}

float3 Renderer::GetEdgeDebugColor(float2 uv)
{
    if (abs(uv.x) < 0.03f || abs(uv.x - 1) < 0.03f || abs(uv.y) < 0.03f || abs(uv.y - 1) < 0.03f)
    {
        return float3(0, 0, 0);
    }
    else { return float3(1); }
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
                    Sample(camera.GetPrimaryRay((float)x + RandomFloat(seed), (float)y + RandomFloat(seed), seed),
                           seed),
                    0);
            float4 pixel = accumulator[x + y * SCRWIDTH] * scale;
            pixel = GammaCorrection(pixel); // Gamma correction, highly decrease frame rate
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
    if (animating) scene.SetTime(anim_time += deltaTime * 0.002f), ClearAccumulator();
    // pixel loop
    Timer t;

    // render tiles using a simple job system
    for (int jobIdx = 0, y = 0; y < SCRHEIGHT / 16; y++)
        for (int x = 0; x < SCRWIDTH / 16; x++)
            jm->AddJob2(tileJob[jobIdx++].Init(this, x, y));
    jm->RunJobs();
    // gather energy received on tiles
    energy = 0;
    for (int tiles = (SCRWIDTH / 16) * (SCRHEIGHT / 16), i = 0; i < tiles; i++)
        energy += tileJob[i].sum;
    // performance report - running average - ms, MRays/s
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
    changed |= ImGui::Checkbox("DOF", &camera.enableDOF);
    changed |= ImGui::SliderFloat("Focal distance", &camera.focalDistance, 1.0f, 10.0f, "%.2f");
    changed |= ImGui::SliderFloat("Aparture size", &camera.aparture, 0.1f, 10.0f, "%.2f");
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
    if (changed) { ClearAccumulator(); }
}