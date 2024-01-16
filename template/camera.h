#pragma once

#include "helper.h"

// default screen resolution
#define SCRWIDTH 1024
#define SCRHEIGHT 640
// #define FULLSCREEN
// #define DOUBLESIZE

namespace Tmpl8
{

class Camera
{
  public:
    Camera()
    {
        // setup a basic view frustum
        camPos = float3(0, 0, -2);
        camTarget = float3(0, 0, -1);
        topLeft = float3(-aspect, 1, 0);
        topRight = float3(aspect, 1, 0);
        bottomLeft = float3(-aspect, -1, 0);
        UpdateCamera();
    }
    Ray GetPrimaryRay(const float x, const float y)
    {
        // calculate pixel position on virtual screen plane
        const float u = (float)x * (1.0f / SCRWIDTH);
        const float v = (float)y * (1.0f / SCRHEIGHT);
        const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
        float3 dir = normalize(P - camPos);

        if (enableDOF)
        {
            float3 fp = camPos + focalDistance * dir;
            float3 randomPoint = UniformRandomPointDisk();
            float3 origin = camPos + randomPoint * aparture;
            return Ray(origin, normalize(fp - origin));
        }
        return Ray(camPos, dir);
    }
    bool HandleInput(const float t)
    {
        if (!WindowHasFocus())
            return false;
        float mSpeed = 0.00025f * t * moveSpeed;
        float tSpeed = 0.00025f * t * turnSpeed;
        float3 ahead = normalize(camTarget - camPos);
        float3 tmpUp(0, 1, 0);
        float3 right = normalize(cross(tmpUp, ahead));
        float3 up = normalize(cross(ahead, right));
        bool changed = false;
        if (IsKeyDown(GLFW_KEY_A))
            camPos -= mSpeed * targetDistance * right, changed = true;
        if (IsKeyDown(GLFW_KEY_D))
            camPos += mSpeed * targetDistance * right, changed = true;
        if (IsKeyDown(GLFW_KEY_W))
            camPos += mSpeed * targetDistance * ahead, changed = true;
        if (IsKeyDown(GLFW_KEY_S))
            camPos -= mSpeed * targetDistance * ahead, changed = true;
        if (IsKeyDown(GLFW_KEY_R))
            camPos += mSpeed * targetDistance * up, changed = true;
        if (IsKeyDown(GLFW_KEY_F))
            camPos -= mSpeed * targetDistance * up, changed = true;
        camTarget = camPos + ahead;
        if (IsKeyDown(GLFW_KEY_UP))
            camTarget -= tSpeed * up, changed = true;
        if (IsKeyDown(GLFW_KEY_DOWN))
            camTarget += tSpeed * up, changed = true;
        if (IsKeyDown(GLFW_KEY_LEFT))
            camTarget -= tSpeed * right, changed = true;
        if (IsKeyDown(GLFW_KEY_RIGHT))
            camTarget += tSpeed * right, changed = true;
        if (!changed)
            return false;
        ahead = normalize(camTarget - camPos);
        up = normalize(cross(ahead, right));
        right = normalize(cross(up, ahead));
        ahead *= targetDistance;
        topLeft = camPos + ahead - aspect * right + up;
        topRight = camPos + ahead + aspect * right + up;
        bottomLeft = camPos + ahead - aspect * right - up;
        return true;
    }
    void UpdateCamera()
    {
        float3 ahead = normalize(camTarget - camPos);
        float3 tmpUp(0, 1, 0);
        float3 right = normalize(cross(tmpUp, ahead));
        float3 up = normalize(cross(ahead, right));
        right = normalize(cross(up, ahead));
        ahead *= targetDistance;
        topLeft = camPos + ahead - aspect * right + up;
        topRight = camPos + ahead + aspect * right + up;
        bottomLeft = camPos + ahead - aspect * right - up;
    }
    void SetCameraState(float3 position, float3 target)
    {
        camPos = position;
        camTarget = target;
        float3 ahead = normalize(camTarget - camPos);
        camTarget = camPos + ahead * targetDistance;
        float3 tmpUp(0, 1, 0);
        float3 right = normalize(cross(tmpUp, ahead));
        float3 up = normalize(cross(ahead, right));
        right = normalize(cross(up, ahead));
        ahead *= targetDistance;
        topLeft = camPos + ahead - aspect * right + up;
        topRight = camPos + ahead + aspect * right + up;
        bottomLeft = camPos + ahead - aspect * right - up;
    }
    float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
    float3 camPos, camTarget;
    float3 topLeft, topRight, bottomLeft;
    bool enableDOF = false;
    float targetDistance = 2;
    float focalDistance = 2;
    float aparture = 1;
    float moveSpeed = 5;
    float turnSpeed = 5;
};

} // namespace Tmpl8