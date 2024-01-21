#include "precomp.h"
#include "mesh.h"

Tmpl8::Mesh::Mesh(const int idx, const std::string &modelPath, mat4 scaleMat)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, modelPath.c_str()))
    {
        throw std::runtime_error(warn + err);
    }

    std::unordered_map<Vertex, uint32_t> uniqueVertices{};

    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;

    for (const auto &shape : shapes)
    {
        for (const auto &index : shape.mesh.indices)
        {
            Vertex vertex{};

            if (index.vertex_index >= 0)
            {
                vertex.position = {attrib.vertices[3 * index.vertex_index + 0],
                                   attrib.vertices[3 * index.vertex_index + 1],
                                   attrib.vertices[3 * index.vertex_index + 2]};
            }

            if (index.normal_index >= 0)
            {
                vertex.normal = {attrib.normals[3 * index.normal_index + 0], attrib.normals[3 * index.normal_index + 1],
                                 attrib.normals[3 * index.normal_index + 2]};
            }

            if (index.texcoord_index >= 0)
            {
                vertex.uv = {attrib.texcoords[2 * index.texcoord_index + 0],
                             attrib.texcoords[2 * index.texcoord_index + 1]};
            }

            if (uniqueVertices.count(vertex) == 0)
            {
                uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
                vertices.push_back(vertex);
            }
            indices.push_back(uniqueVertices[vertex]);
        }
    }

    meshIdx = idx;

    triCount = indices.size();

    triangles = new Tri[triCount];
    triangleExs = new TriEx[triCount];

    for (int i = 0; i < triCount; i += 3)
    {
        Tri tri(TransformPosition(vertices[indices[i]].position, scaleMat),
                TransformPosition(vertices[indices[i + 1]].position, scaleMat),
                TransformPosition(vertices[indices[i + 2]].position, scaleMat), float3(0));
        TriEx triEx(vertices[indices[i]].normal, vertices[indices[i + 1]].normal, vertices[indices[i + 2]].normal,
                    vertices[indices[i]].uv, vertices[indices[i + 1]].uv, vertices[indices[i + 2]].uv);
        tri.centroid = (tri.vertex0 + tri.vertex1 + tri.vertex2) * 0.3333f;
        triangles[i] = tri;
        triangleExs[i] = triEx;
    }
}
