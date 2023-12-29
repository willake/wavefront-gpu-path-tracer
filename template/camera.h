#pragma once

// default screen resolution
#define SCRWIDTH	1024
#define SCRHEIGHT	640
// #define FULLSCREEN
// #define DOUBLESIZE

namespace Tmpl8 {

class Camera
{
public:
	Camera()
	{
		// setup a basic view frustum
		camPos = float3( 0, 0, -2 );
		camTarget = float3( 0, 0, -1 );
		topLeft = float3( -aspect, 1, 0 );
		topRight = float3( aspect, 1, 0 );
		bottomLeft = float3( -aspect, -1, 0 );
	}
	Ray GetPrimaryRay( const float x, const float y )
	{
		// calculate pixel position on virtual screen plane
		const float u = (float)x * (1.0f / SCRWIDTH);
		const float v = (float)y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		return Ray( camPos, normalize( P - camPos ) );
	}
	bool HandleInput( const float t )
	{
		if (!WindowHasFocus()) return false;
		float mSpeed = 0.00025f * t * moveSpeed;
		float tSpeed = 0.00025f * t * turnSpeed;
		float3 ahead = normalize( camTarget - camPos );
		float3 tmpUp( 0, 1, 0 );
		float3 right = normalize( cross( tmpUp, ahead ) );
		float3 up = normalize( cross( ahead, right ) );
		bool changed = false;
		if (IsKeyDown( GLFW_KEY_A )) camPos -= mSpeed * 2 * right, changed = true;
		if (IsKeyDown( GLFW_KEY_D )) camPos += mSpeed * 2 * right, changed = true;
		if (IsKeyDown( GLFW_KEY_W )) camPos += mSpeed * 2 * ahead, changed = true;
		if (IsKeyDown( GLFW_KEY_S )) camPos -= mSpeed * 2 * ahead, changed = true;
		if (IsKeyDown( GLFW_KEY_R )) camPos += mSpeed * 2 * up, changed = true;
		if (IsKeyDown( GLFW_KEY_F )) camPos -= mSpeed * 2 * up, changed = true;
		camTarget = camPos + ahead;
		if (IsKeyDown( GLFW_KEY_UP )) camTarget -= tSpeed * up, changed = true;
		if (IsKeyDown( GLFW_KEY_DOWN )) camTarget += tSpeed * up, changed = true;
		if (IsKeyDown( GLFW_KEY_LEFT )) camTarget -= tSpeed * right, changed = true;
		if (IsKeyDown( GLFW_KEY_RIGHT )) camTarget += tSpeed * right, changed = true;
		if (!changed) return false;
		ahead = normalize( camTarget - camPos );
		up = normalize( cross( ahead, right ) );
		right = normalize( cross( up, ahead ) );
		topLeft = camPos + 2 * ahead - aspect * right + up;
		topRight = camPos + 2 * ahead + aspect * right + up;
		bottomLeft = camPos + 2 * ahead - aspect * right - up;
		return true;
	}
	void SetCameraState(float3 position, float3 target)
	{
		camPos = position;
		camTarget = target;
		float3 ahead = normalize(camTarget - camPos);
		float3 tmpUp(0, 1, 0);
		float3 right = normalize(cross(tmpUp, ahead));
		float3 up = normalize(cross(ahead, right));
		right = normalize(cross(up, ahead));
		topLeft = camPos + 2 * ahead - aspect * right + up;
		topRight = camPos + 2 * ahead + aspect * right + up;
		bottomLeft = camPos + 2 * ahead - aspect * right - up;
	}
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float3 camPos, camTarget;
	float3 topLeft, topRight, bottomLeft;
	float moveSpeed = 5;
	float turnSpeed = 5;
};

}