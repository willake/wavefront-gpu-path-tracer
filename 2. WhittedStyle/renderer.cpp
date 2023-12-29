#include "precomp.h"
#include "helper.h"
#include "renderer.h"

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
	m_buildTime = scene.GetBuildTime();
	m_triangleCount = scene.GetTriangleCount();
	m_maxTreeDepth = scene.GetMaxTreeDepth();
	kernelGeneratePrimaryRays = Kernel("../assets/cl/kernels.cl", "generatePrimaryRays");
	kernelExtend = Kernel("../assets/cl/kernels.cl", "extend");
	kernelShade = Kernel("../assets/cl/kernels.cl", "shade");
	kernelConnect = Kernel("../assets/cl/kernels.cl", "connect");
	kernelFinalize = Kernel("../assets/cl/kernels.cl", "finalize");
	kernelGeneratePrimaryRays.SetArguments(rayBuffer);
	/*rayBuffer = Buffer(SCRWIDTH * SCRHEIGHT * sizeof(Ray));*/
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace(Ray& ray, int depth)
{
	if (depth > depthLimit) return float3(0);
	scene.FindNearest(ray);
	//if (ray.objIdx == -1) return float3(0);
	if (ray.objIdx == -1) return scene.GetSkyColor(ray); // or a fancy sky color
	float3 I = ray.O + ray.t * ray.D;
	HitInfo hitInfo = scene.GetHitInfo(ray, I);
	float3 N = hitInfo.normal;
	float2 uv = hitInfo.uv;
	Material* material = hitInfo.material;
	float3 albedo = material->isAlbedoOverridden ? scene.GetAlbedo(ray.objIdx, I) : material->GetAlbedo(uv);

	/* visualize edges */ // return GetEdgeDebugColor(ray.barycentric);
	/* visualize normal */ // return N; // return (N + 1) * 0.5f;
	/* visualize distance */ // return 0.1f * float3( ray.t, ray.t, ray.t );
	/* visualize albedo */ // return albedo;
	if (m_inspectTraversal) return GetTraverseCountColor(ray.traversed, m_peakTraversal);
	if (m_inspectIntersectionTest) return GetTraverseCountColor(ray.tested, m_peakTests);

	if (material->isLight) return scene.GetLightColor();

	float3 out_radiance(0);
	float reflectivity = material->reflectivity;
	float refractivity = material->refractivity;
	float diffuseness = 1 - (reflectivity + refractivity);

	if (reflectivity > 0.0f)
	{
		float3 R = reflect(ray.D, N);
		Ray r(I + R * EPSILON, R);
		out_radiance += reflectivity * albedo * Trace(r, depth + 1);
	}
	else if (refractivity > 0.0f)
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
			out_radiance += albedo * (1 - Fr) * Trace(t, depth + 1);
		}
		out_radiance += albedo * Fr * Trace(r, depth + 1);
	}

	if (diffuseness > 0)
	{
		float3 irradiance = DirectIllumination(I, N);
		float3 ambient = float3(0.3f, 0.3f, 0.3f);
		float3 brdf = albedo * INVPI;
		out_radiance += diffuseness * brdf * (irradiance + ambient);
	}
	float3 medium_scale(1);
	if (ray.inside)
	{
		float3 absorption = material->absorption;
		medium_scale.x = expf(absorption.x * -ray.t);
		medium_scale.y = expf(absorption.y * -ray.t);
		medium_scale.z = expf(absorption.z * -ray.t);
	}

	return medium_scale * out_radiance;
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

float3 Renderer::DirectIllumination(float3 I, float3 N)
{
	// sum irradiance from light sources
	float3 irradiance(0);
	// query the (only) scene light
	float3 pointOnLight = scene.GetLightPos();
	float3 L = pointOnLight - I;
	float distance = length(L);
	L *= 1 / distance;
	float ndotl = dot(N, L);
	if (ndotl < EPSILON) /* we don't face the light */ return 0;
	// cast a shadow ray
	Ray s(I + L * EPSILON, L, distance - 2 * EPSILON);
	if (!scene.IsOccluded(s))
	{
		// light is visible; calculate irradiance (= projected radiance)
		float attenuation = 1 / (distance * distance);
		float3 in_radiance = scene.GetLightColor() * attenuation;
		irradiance = in_radiance * dot(N, L);
	}
	return irradiance;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
	// animation
	if (animating) scene.SetTime(anim_time += deltaTime * 0.002f);
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			Ray primaryRay = camera.GetPrimaryRay((float)x, (float)y);
			float4 pixel = float4(Trace(primaryRay, 0), 0);

			// for metrics
			if (primaryRay.traversed > 0) m_rayHitCount++;
			if (primaryRay.traversed > m_peakTraversal) m_peakTraversal = primaryRay.traversed;
			if (primaryRay.tested > m_peakTests) m_peakTests = primaryRay.tested;
			m_totalTraversal += primaryRay.traversed;
			m_totalTests += primaryRay.tested;
			// translate accumulator contents to rgb32 pixels
			screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
			accumulator[x + y * SCRWIDTH] = pixel;
		}
	}
	// performance report - running average - ms, MRays/s
	/*static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf("%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000);*/
	if (m_rayHitCount > 0)
	{
		m_averageTraversal = m_totalTraversal / m_rayHitCount;
		m_averageTests = m_totalTests / m_rayHitCount;
	}
	m_avg = (1 - m_alpha) * m_avg + m_alpha * t.elapsed() * 1000;
	if (m_alpha > 0.05f) m_alpha *= 0.5f;
	m_fps = 1000.0f / m_avg, m_rps = (SCRWIDTH * SCRHEIGHT) / m_avg;
	system("cls");
	printf("Total Traversal: %.0f\n", m_totalTraversal);
	printf("Average Traversal: %.2f\n", m_averageTraversal);
	printf("Peak Traversal: %.2f\n", m_peakTraversal);
	printf("Total Tests: %.0f\n", m_totalTests);
	printf("Average Tests: %.2f\n", m_averageTests);
	printf("Peak Tests: %.2f\n", m_peakTests);
	// handle user input
	if (camera.HandleInput(deltaTime))
	{
		m_averageTraversal = 0;
		m_averageTests = 0;
		m_peakTraversal = 0;
		m_peakTests = 0;
	}
	m_totalTraversal = 0;
	m_totalTests = 0;
	m_rayHitCount = 0;
}

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
	// animation toggle
	ImGui::Checkbox("Animate scene", &animating);
	ImGui::Checkbox("Inspect traversal", &m_inspectTraversal);
	ImGui::Checkbox("Inspect intersection", &m_inspectIntersectionTest);
	ImGui::SliderFloat("Camera move speed", &camera.moveSpeed, 1.0f, 10.0f, "%.2f");
	ImGui::SliderFloat("Camera turn speed", &camera.turnSpeed, 1.0f, 10.0f, "%.2f");
	// camera position field
	ImGuiFloat3("Position", m_camPositionToSet);
	ImGuiFloat3("Rotation", m_camTargetToSet);
	if (ImGui::Button("Set Camera")) {
		// Button was clicked, perform action (e.g., reset values)
		camera.SetCameraState(m_camPositionToSet, m_camTargetToSet);
		m_averageTraversal = 0;
		m_averageTests = 0;
		m_peakTraversal = 0;
		m_peakTests = 0;
	}
	// ray query on mouse
	Ray r = camera.GetPrimaryRay((float)mousePos.x, (float)mousePos.y);
	scene.FindNearest(r);
	ImGui::Text("Object id: %i", r.objIdx);
	ImGui::Text("Triangle count: %i", m_triangleCount);
	ImGui::Text("Build time: %lld", m_buildTime.count());
	ImGui::Text("Max tree depth: %d", m_maxTreeDepth);
	ImGui::Text("Frame: %5.2f ms (%.1ffps)", m_avg, m_fps);
	//ImGui::Text("FPS: %.1ffps", m_fps);
	ImGui::Text("RPS: %.1f Mrays/s", m_rps);
	ImGui::Text("Camera pos: (%.2f, %.2f, %.2f)", camera.camPos.x, camera.camPos.y, camera.camPos.z);
	ImGui::Text("Camera target: (%.2f, %.2f, %.2f)", camera.camTarget.x, camera.camTarget.y, camera.camTarget.z);
}