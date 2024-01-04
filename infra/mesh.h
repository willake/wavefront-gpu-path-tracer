#pragma once

namespace Tmpl8
{
	struct MeshInstance
	{
		int meshIdx = -1;
		int startIdx = -1;
	};

	class Mesh
	{
	public:
		Mesh(const int idx, const std::string& modelPath, mat4 transform);
	public:
		int meshIdx = -1;
		std::vector<Tri> triangles;
		std::vector<TriEx> triangleExs;
	};
}
