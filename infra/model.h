#pragma once

namespace Tmpl8
{
	class Model
	{
	public:
		Model() {}
		//Model(const int idx, const std::string& path, mat4 transform);
		Model(const int idx, const std::string& modelPath, mat4 transform);
		void AppendTriangles(std::vector<Tri>& triangles, std::vector<TriEx>& triangleExs);
	public:
		int objIdx = -1;
		int matIdx = -1;
		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		mat4 T, invT;
	};
}