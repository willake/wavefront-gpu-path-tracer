#pragma once
namespace Tmpl8
{
// contain all buffers of a scene
class SceneBuffer
{
  public:
    void CreateBuffers();
    void CopyToDevice();
    uint objCount = 0;
    uint materialCount = 0;
    uint meshCount = 0;
    uint lightCount = 0;
    uint totalTriangleCount = 0;
    uint totalBVHNodeCount = 0;
    uint totalPixelCount = 0;

    uint *skydomePixels;
    GPUTexture skydomeTexture;
    Buffer *skydomeBuffer;
    uint *floorPixels;
    GPUTexture floorTexture;
    Buffer *floorBuffer;
    Tri *triangles;
    Buffer *triBuffer;
    TriEx *triangleExs;
    Buffer *triExBuffer;
    uint *triangleIndices;
    Buffer *triIdxBuffer;
    MeshInstance *meshInstances;
    Buffer *meshInsBuffer;
    BVHNode *bvhNodes;
    Buffer *bvhNodeBuffer;
    GPUBVH *gpubvhs;
    Buffer *bvhBuffer;
    GPUBLAS *gpublases;
    Buffer *blasBuffer;
    TLASNode *tlasNodes;
    Buffer *tlasNodeBuffer;
    uint *texturePixels;
    Buffer *texturePixelBuffer;
    GPUTexture *gputextures;
    Buffer *textureBuffer;
    GPUMaterial *gpuMats;
    Buffer *materialBuffer;
    Light *lights;
    Buffer *lightBuffer;
};
} // namespace Tmpl8