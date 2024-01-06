#include "precomp.h"
#include "tlas_file_scene.h"

float3 TLASFileScene::RandomPointOnLight(const float r0, const float r1, uint &lightIdx) const
{
    //// select a random light and use that
    lightIdx = (uint)(r0 * lightCount);
    const Quad &q = lightQuads[lightIdx];
    // renormalize r0 for reuse
    float stratum = lightIdx * 0.25f;
    float r2 = (r0 - stratum) / (1 - stratum);
    // get a random position on the selected quad
    const float size = q.size;
    float3 corner1 = TransformPosition(float3(-size, 0, -size), q.T);
    float3 corner2 = TransformPosition(float3(size, 0, -size), q.T);
    float3 corner3 = TransformPosition(float3(-size, 0, size), q.T);
    return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
}

TLASFileScene::TLASFileScene(const string &filePath)
{
    errorMaterial.albedo = float3(255, 192, 203) / 255.f;

    SceneData sceneData = LoadSceneFile(filePath);

    primitiveMaterials[0].isLight = true;
    // materials[1].isAlbedoOverridden = true;
    // primitiveMaterials[1].textureDiffuse = std::make_unique<Texture>(sceneData.planeTextureLocation);
    primitiveMaterials[1].textureDiffuse = Texture(sceneData.planeTextureLocation);
    objIdUsed = 2;

    floor = Plane(1, float3(0, 1, 0), 1, primitiveMaterials[1].textureDiffuse.width / 100);

    sceneName = sceneData.name;
    skydome = Texture(sceneData.skydomeLocation);

    // create lights
    lightCount = sceneData.lights.size();

    lightQuads = new Quad[lightCount];
    lights = new Light[lightCount];

    for (int i = 0; i < lightCount; i++)
    {
        LightData &lightData = sceneData.lights[i];
        lightQuads[i] = Quad(0, lightData.size);
        mat4 T = mat4::Translate(lightData.position);
        // *mat4::RotateZ(sinf(animTime * 0.6f) * 0.1f);
        lightQuads[i].T = T, lightQuads[i].invT = T.FastInvertedTransformNoScale();
        lights[i] = Light();
        lights[i].position = lightData.position;
        lights[i].normal = lightQuads[i].GetNormal(float3(0));
        lights[i].area = lightData.size * lightData.size;
    }

    objCount = sceneData.objects.size();

    materialCount = sceneData.materials.size();

    materials = new Material[materialCount];
    gpuMats = new GPUMaterial[materialCount];

    meshCount = sceneData.meshes.size();

    // Setup materials
    for (int i = 0; i < materialCount; i++)
    {
        materials[i] = Material();
        materials[i].reflectivity = sceneData.materials[i].reflectivity;
        materials[i].refractivity = sceneData.materials[i].refractivity;
        materials[i].absorption = sceneData.materials[i].absorption;

        gpuMats[i] = GPUMaterial();
        gpuMats[i].reflectivity = sceneData.materials[i].reflectivity;
        gpuMats[i].refractivity = sceneData.materials[i].refractivity;
        gpuMats[i].absorption = sceneData.materials[i].absorption;

        if (!sceneData.materials[i].textureLocation.empty())
            materials[i].textureDiffuse = Texture(sceneData.materials[i].textureLocation);
        // materials[i].textureDiffuse = std::make_unique<Texture>(sceneData.materials[i].textureLocation);
    }

    // Setup texturePixels

    gputextures = new GPUTexture[materialCount];
    totalPixelCount = 0;
    for (int i = 0; i < materialCount; i++)
    {
        Texture &texture = materials[i].textureDiffuse;
        gputextures[i].width = texture.width;
        gputextures[i].height = texture.height;
        gputextures[i].startIdx = totalPixelCount;
        totalPixelCount += texture.width * texture.height;
    }

    texturePixels = new uint[totalPixelCount];

    for (int i = 0; i < materialCount; i++)
    {
        Texture &texture = materials[i].textureDiffuse;
        uint pixelCount = texture.width * texture.height;

        for (int pI = 0; pI < pixelCount; pI++)
        {
            texturePixels[gputextures[i].startIdx + pI] = texture.pixels[pI];
        }
    }
    // Setup meshes

    meshes.resize(meshCount);
    meshInstances = new MeshInstance[meshCount];

    // prepare for the huge triangle array
    totalTriangleCount = 0;

    for (int i = 0; i < meshCount; i++)
    {
        MeshData &meshData = sceneData.meshes[i];
        mat4 S = mat4::Scale(meshData.scale);
        meshes[i] = Mesh(i, meshData.modelLocation, S);
        totalTriangleCount += meshes[i].triCount;
    }

    triangles = new Tri[totalTriangleCount];
    triangleExs = new TriEx[totalTriangleCount];

    int tmpTriangleIndx = 0;
    for (int i = 0; i < meshCount; i++)
    {
        Mesh &mesh = meshes[i];

        for (int tI = 0; tI < mesh.triCount; tI++)
        {
            triangles[tmpTriangleIndx + tI] = mesh.triangles[tI];
            triangleExs[tmpTriangleIndx + tI] = mesh.triangleExs[tI];
        }

        meshInstances[i].meshIdx = mesh.meshIdx;
        meshInstances[i].triStartIdx = tmpTriangleIndx;
        meshInstances[i].triCount = mesh.triCount;

        tmpTriangleIndx += mesh.triCount;
    }

    // Setup BVHs

    totalBVHNodeCount = 0;
    // prepare for bvh nodes
    bvhs = new BVH[meshCount];
    gpubvhs = new GPUBVH[meshCount];

    for (int i = 0; i < meshCount; i++)
    {
        bvhs[i] = BVH(meshInstances[i], triangles, triangleExs);
        totalBVHNodeCount += bvhs[i].triangleCount * 2 - 1;
    }

    bvhNodes = new BVHNode[totalBVHNodeCount];
    triangleIndices = new uint[totalTriangleCount];

    int tmpBVHNodeIdx = 0;
    for (int i = 0; i < meshCount; i++)
    {
        // put all nodes
        int nodeCount = bvhs[i].triangleCount * 2 - 1;
        for (int bI = 0; bI < nodeCount; bI++)
        {
            bvhNodes[tmpBVHNodeIdx + bI] = bvhs[i].bvhNodes[bI];
        }

        MeshInstance &meshIns = meshInstances[i];
        // put all indices
        for (int tI = 0; tI < bvhs[i].triangleCount; tI++)
        {
            triangleIndices[meshIns.triStartIdx + tI] = bvhs[i].triangleIndices[tI];
        }

        gpubvhs[i] = GPUBVH(i, tmpBVHNodeIdx, nodeCount);
        tmpBVHNodeIdx += gpubvhs[i].nodeCount;
    }

    // Setup BLASes
    blases = new BLAS[objCount];
    gpublases = new GPUBLAS[objCount];

    for (int i = 0; i < objCount; i++)
    {
        ObjectData &objectData = sceneData.objects[i];
        mat4 T = mat4::Translate(objectData.position) * mat4::RotateX(objectData.rotation.x * Deg2Red) *
                 mat4::RotateY(objectData.rotation.y * Deg2Red) * mat4::RotateZ(objectData.rotation.z * Deg2Red);
        blases[i] = BLAS(objIdUsed, &bvhs[objectData.meshIdx], objectData.materialIdx, T);
        gpublases[i] =
            GPUBLAS(objIdUsed, objectData.meshIdx, objectData.materialIdx, T, blases[i].aabbMin, blases[i].aabbMax);
        objIdUsed++;
    }

    // setup tlas
    tlas = TLAS(blases, objCount);

    PrepareBuffers();
    SetTime(0);
}

SceneData TLASFileScene::LoadSceneFile(const string &filePath)
{
    SceneData sceneData;
    // Read the XML file into a string
    std::ifstream file(filePath.c_str());
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Parse the XML content using RapidXML
    rapidxml::xml_document<> doc;
    doc.parse<0>(&content[0]);

    // Get the root node
    rapidxml::xml_node<> *root = doc.first_node("scene");

    // Extract scene information
    sceneData.name = root->first_node("scene_name")->value();
    sceneData.planeTextureLocation = root->first_node("plane_texture_location")->value();
    sceneData.skydomeLocation = root->first_node("skydome_location")->value();

    // Extract light information
    for (rapidxml::xml_node<> *lightNode = root->first_node("lights")->first_node("light"); lightNode;
         lightNode = lightNode->next_sibling())
    {
        LightData light;
        for (rapidxml::xml_node<> *posNode = lightNode->first_node("position")->first_node(); posNode;
             posNode = posNode->next_sibling())
        {
            int index = posNode->name()[0] - 'x'; // 'x', 'y', 'z' map to 0, 1, 2
            light.position[index] = std::stof(posNode->value());
        }
        light.size = std::stoi(lightNode->first_node("size")->value());
        sceneData.lights.push_back(light);
    }

    // Extract object information
    for (rapidxml::xml_node<> *objNode = root->first_node("objects")->first_node("object"); objNode;
         objNode = objNode->next_sibling())
    {
        ObjectData obj;
        obj.meshIdx = std::stoi(objNode->first_node("mesh_idx")->value());
        obj.materialIdx = std::stoi(objNode->first_node("material_idx")->value());

        for (rapidxml::xml_node<> *posNode = objNode->first_node("position")->first_node(); posNode;
             posNode = posNode->next_sibling())
        {
            int index = posNode->name()[0] - 'x'; // 'x', 'y', 'z' map to 0, 1, 2
            obj.position[index] = std::stof(posNode->value());
        }

        for (rapidxml::xml_node<> *rotNode = objNode->first_node("rotation")->first_node(); rotNode;
             rotNode = rotNode->next_sibling())
        {
            int index = rotNode->name()[0] - 'x'; // 'x', 'y', 'z' map to 0, 1, 2
            obj.rotation[index] = std::stof(rotNode->value());
        }

        sceneData.objects.push_back(obj);
    }

    // Extract mesh information
    for (rapidxml::xml_node<> *meshNode = root->first_node("meshes")->first_node("mesh"); meshNode;
         meshNode = meshNode->next_sibling())
    {
        MeshData mesh;
        mesh.modelLocation = meshNode->first_node("model_location")->value();
        for (rapidxml::xml_node<> *scaleNode = meshNode->first_node("scale")->first_node(); scaleNode;
             scaleNode = scaleNode->next_sibling())
        {
            int index = scaleNode->name()[0] - 'x'; // 'x', 'y', 'z' map to 0, 1, 2
            mesh.scale[index] = std::stof(scaleNode->value());
        }

        sceneData.meshes.push_back(mesh);
    }

    // Extract material information
    for (rapidxml::xml_node<> *matNode = root->first_node("materials")->first_node("material"); matNode;
         matNode = matNode->next_sibling())
    {
        MaterialData material;
        material.reflectivity = std::stof(matNode->first_node("reflectivity")->value());
        material.refractivity = std::stof(matNode->first_node("refractivity")->value());

        for (rapidxml::xml_node<> *absNode = matNode->first_node("absorption")->first_node(); absNode;
             absNode = absNode->next_sibling())
        {
            int index = absNode->name()[0] - 'x'; // 'x', 'y', 'z' map to 0, 1, 2
            material.absorption[index] = std::stof(absNode->value());
        }

        material.textureLocation = matNode->first_node("texture_location")->value();

        sceneData.materials.push_back(material);
    }

    return sceneData;
}

void TLASFileScene::PrepareBuffers()
{
    // skydome
    skydomeBuffer = new Buffer(skydome.width * skydome.height * sizeof(uint), skydome.pixels);
    skydomeBuffer->CopyToDevice(true);
    // floor
    Texture &floorTex = primitiveMaterials[1].textureDiffuse;
    floorBuffer = new Buffer(floorTex.width * floorTex.height * sizeof(uint), floorTex.pixels);
    floorBuffer->CopyToDevice(true);
    // materials
    materialBuffer = new Buffer(materialCount * sizeof(GPUMaterial), gpuMats);
    materialBuffer->CopyToDevice(true);
    texturePixelBuffer = new Buffer(totalPixelCount * sizeof(uint), texturePixels);
    texturePixelBuffer->CopyToDevice(true);
    textureBuffer = new Buffer(materialCount * sizeof(GPUTexture), gputextures);
    textureBuffer->CopyToDevice(true);
    // scene data
    triBuffer = new Buffer(totalTriangleCount * sizeof(Tri), triangles);
    triBuffer->CopyToDevice(true);
    triExBuffer = new Buffer(totalTriangleCount * sizeof(TriEx), triangleExs);
    triExBuffer->CopyToDevice(true);
    triIdxBuffer = new Buffer(totalTriangleCount * sizeof(uint), triangleIndices);
    triIdxBuffer->CopyToDevice(true);
    meshInsBuffer = new Buffer(meshCount * sizeof(MeshInstance), meshInstances);
    meshInsBuffer->CopyToDevice(true);
    bvhNodeBuffer = new Buffer(totalBVHNodeCount * sizeof(BVHNode), bvhNodes);
    bvhNodeBuffer->CopyToDevice(true);
    bvhBuffer = new Buffer(meshCount * sizeof(GPUBVH), gpubvhs);
    bvhBuffer->CopyToDevice(true);
    blasBuffer = new Buffer(objCount * sizeof(GPUBLAS), gpublases);
    blasBuffer->CopyToDevice(true);
    tlasNodeBuffer = new Buffer(objCount * 2 * sizeof(TLASNode), tlas.tlasNode);
    tlasNodeBuffer->CopyToDevice(true);
}

void TLASFileScene::SetTime(float t)
{
    animTime = t;
}

float3 TLASFileScene::GetSkyColor(const Ray &ray) const
{
    // Convert ray direction to texture coordinates, assuming a spherical skydome
    float phi = atan2(-ray.D.z, ray.D.x) + PI;
    float theta = acos(-ray.D.y);
    float u = phi * INV2PI;
    float v = theta * INVPI;

    //// Sample the HDR skydome texture
    float3 color = skydome.Sample(u, v);

    return color;
}

Light Tmpl8::TLASFileScene::GetLight(int idx)
{
    if (idx >= lightCount)
        return lights[0];
    return lights[idx];
}

float3 TLASFileScene::GetLightPos() const
{
    // light point position is the middle of the swinging quad
    float3 corner1 = TransformPosition(float3(-0.5f, 0, -0.5f), lightQuads[0].T);
    float3 corner2 = TransformPosition(float3(0.5f, 0, 0.5f), lightQuads[0].T);
    return (corner1 + corner2) * 0.5f - float3(0, 0.01f, 0);
}

float3 TLASFileScene::GetLightColor() const
{
    return float3(24, 24, 22);
}

float3 TLASFileScene::RandomPointOnLight(uint &seed, uint &lightIdx) const
{
    return RandomPointOnLight(RandomFloat(seed), RandomFloat(seed), lightIdx);
}

void TLASFileScene::FindNearest(Ray &ray)
{
    for (int i = 0; i < lightCount; i++)
    {
        lightQuads[i].Intersect(ray);
    }
    floor.Intersect(ray);
    tlas.Intersect(ray);
}

bool TLASFileScene::IsOccluded(const Ray &ray)
{
    for (int i = 0; i < lightCount; i++)
    {
        if (lightQuads[i].IsOccluded(ray))
            return true;
    }
    Ray shadow = Ray(ray);
    shadow.t = 1e34f;
    tlas.Intersect(shadow);
    if (shadow.objIdx > -1)
        return true;
    // skip planes
    return false;
}

HitInfo TLASFileScene::GetHitInfo(const Ray &ray, const float3 I)
{
    HitInfo hitInfo = HitInfo(float3(0), float2(0), &errorMaterial);
    switch (ray.objIdx)
    {
    case 0:
        hitInfo.normal = lightQuads[0].GetNormal(I);
        hitInfo.uv = float2(0);
        hitInfo.material = &primitiveMaterials[0];
        break;
    case 1:
        hitInfo.normal = floor.GetNormal(I);
        hitInfo.uv = floor.GetUV(I);
        hitInfo.material = &primitiveMaterials[1];
        break;
    default:
        BLAS &blas = tlas.blases[ray.objIdx - 2];
        hitInfo.normal = blas.GetNormal(ray.triIdx, ray.barycentric);
        hitInfo.uv = blas.GetUV(ray.triIdx, ray.barycentric);
        hitInfo.material = &materials[blas.matIdx];
        break;
    }

    if (dot(hitInfo.normal, ray.D) > 0)
        hitInfo.normal = -hitInfo.normal;

    return hitInfo;
}

float3 TLASFileScene::GetAlbedo(int objIdx, float3 I) const
{
    if (objIdx == 1)
        return floor.GetAlbedo(I);
    return float3(0);
}

int TLASFileScene::GetTriangleCount() const
{
    return totalTriangleCount;
}

std::chrono::microseconds TLASFileScene::GetBuildTime() const
{
    std::chrono::microseconds time(0);
    for (int i = 0; i < meshCount; i++)
    {
        time += bvhs[i].buildTime;
    }
    time += tlas.buildTime;
    return time;
}

uint TLASFileScene::GetMaxTreeDepth() const
{
    uint maxDepth = 0;
    for (int i = 0; i < meshCount; i++)
    {
        if (bvhs[i].maxDepth > maxDepth)
            maxDepth = bvhs[i].maxDepth;
    }
    return maxDepth;
}