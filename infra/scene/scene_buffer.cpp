#include "precomp.h"
#include "scene_buffer.h"

void SceneBuffer::CreateBuffers()
{
    // skydome
    skydomeBuffer = new Buffer(sceneProperty.skydomeTexture.width * sceneProperty.skydomeTexture.height * sizeof(uint),
                               skydomePixels);
    // floor
    floorBuffer =
        new Buffer(sceneProperty.floorTexture.width * sceneProperty.floorTexture.height * sizeof(uint), floorPixels);
    // materials
    materialBuffer = new Buffer(materialCount * sizeof(GPUMaterial), gpuMats);
    texturePixelBuffer = new Buffer(totalPixelCount * sizeof(uint), texturePixels);
    textureBuffer = new Buffer(materialCount * sizeof(GPUTexture), gputextures);
    // scene data
    triBuffer = new Buffer(totalTriangleCount * sizeof(Tri), triangles);
    triExBuffer = new Buffer(totalTriangleCount * sizeof(TriEx), triangleExs);
    triIdxBuffer = new Buffer(totalTriangleCount * sizeof(uint), triangleIndices);
    meshInsBuffer = new Buffer(meshCount * sizeof(MeshInstance), meshInstances);
    bvhNodeBuffer = new Buffer(totalBVHNodeCount * sizeof(BVHNode), bvhNodes);
    bvhBuffer = new Buffer(meshCount * sizeof(GPUBVH), gpubvhs);
    blasBuffer = new Buffer(objCount * sizeof(GPUBLAS), gpublases);
    tlasNodeBuffer = new Buffer(objCount * 2 * sizeof(TLASNode), tlasNodes);
    lightBuffer = new Buffer(lightCount * sizeof(Light), lights);
    scenePropertyBuffer = new Buffer(sizeof(SceneProperty), &sceneProperty);
}

void SceneBuffer::CopyToDevice()
{
    skydomeBuffer->CopyToDevice(true);
    floorBuffer->CopyToDevice(true);
    materialBuffer->CopyToDevice(true);
    texturePixelBuffer->CopyToDevice(true);
    textureBuffer->CopyToDevice(true);
    triBuffer->CopyToDevice(true);
    triExBuffer->CopyToDevice(true);
    triIdxBuffer->CopyToDevice(true);
    meshInsBuffer->CopyToDevice(true);
    bvhNodeBuffer->CopyToDevice(true);
    bvhBuffer->CopyToDevice(true);
    blasBuffer->CopyToDevice(true);
    tlasNodeBuffer->CopyToDevice(true);
    lightBuffer->CopyToDevice(true);
    scenePropertyBuffer->CopyToDevice(true);
}