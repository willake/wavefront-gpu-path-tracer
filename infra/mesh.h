#pragma once

namespace Tmpl8
{
	struct MeshInstance
	{
		int meshIdx = -1; // 4 bytes
		int triStartIdx = -1; // 4 bytes
		int triCount = 0; // 4 bytes 
		// 12 bytes in total
	};

	class Mesh
	{
	public:
		Mesh() {}
		Mesh(const int idx, const std::string& modelPath, mat4 transform);
	public:
		int meshIdx = -1;
		int triCount = 0;
		Tri* triangles = nullptr;
		TriEx* triangleExs = nullptr;
	};
}
