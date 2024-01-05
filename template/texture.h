#pragma once

// #define STB_IMAGE_IMPLEMENTATION
// #include "stb_image.h"

// some of codes copy from surface.h

#include <string>

namespace Tmpl8
{
__declspec(align(16)) struct GPUTexture
{
    GPUTexture(){};
    GPUTexture(uint w, uint h) : width(w), height(h){};
    uint width = 0;    // 4 bytes
    uint height = 0;   // 4 bytes
    uint startIdx = 0; // 4 bytes
    uint dummy = 0;    // 4 bytes
                       // 16 bytes in total
};
class Texture
{
  private:
    void LoadFromFile(const std::string &file)
    {
        int n;
        unsigned char *data = stbi_load(file.c_str(), &width, &height, &n, 0);
        if (data)
        {
            // pixels = (uint*)MALLOC64(width * height * sizeof(uint));
            pixels = new uint[width * height];
            // pixels.resize(width * height);
            ownBuffer = true; // needs to be deleted in destructor
            const int s = width * height;
            if (n == 1) // greyscale
            {
                for (int i = 0; i < s; i++)
                {
                    const unsigned char p = data[i];
                    pixels[i] = p + (p << 8) + (p << 16);
                }
            }
            else
            {
                for (int i = 0; i < s; i++)
                {
                    pixels[i] = (data[i * n + 0] << 16) + (data[i * n + 1] << 8) + data[i * n + 2];
                }
            }
        }
        stbi_image_free(data);
    }

  public:
    Texture() = default;
    Texture(const std::string &file) : pixels(0), width(0), height(0), ownBuffer(false)
    {
        FILE *f = fopen(file.c_str(), "rb");
        if (!f)
            FatalError("File not found: %s", file);
        fclose(f);
        LoadFromFile(file);
    }
    // Texture(const Texture& texture)
    //{
    //	pixels = texture.pixels;
    //	//pixels.insert(pixels.begin(), texture.pixels.begin(), texture.pixels.end());
    //	int width = texture.width, height = texture.height;
    //	bool ownBuffer = texture.ownBuffer;
    // }
    ~Texture()
    {
        // if (ownBuffer) FREE64(pixels); // free only if we allocated the buffer ourselves
    }

    float3 Sample(float u, float v) const
    {
        if (!pixels)
        // if (pixels.size() == 0)
        {
            // throw exception (texture not loaded)
            return float3(0);
        }
        // Clamp texture coordinates to [0, 1]
        // u = fmod(u, 1.0f);
        // v = fmod(v, 1.0f);

        u = clamp(u, 0.0f, 1.0f);
        v = 1 - clamp(v, 0.0f, 1.0f);

        // Calculate pixel coordinates
        int x = static_cast<int>(u * width);
        int y = static_cast<int>(v * height);

        // Ensure coordinates are within bounds
        x = clamp(x, 0, width - 1);
        y = clamp(y, 0, height - 1);

        // Calculate index in the image data array
        int index = x + y * width;

        uint pixel = pixels[index];

        // Sample color from the texture
        float rgbScale = 1 / 255.0f;
        float r = ((pixel >> 16) & 0xFF) * rgbScale;
        float g = ((pixel >> 8) & 0xFF) * rgbScale;
        float b = (pixel & 0xFF) * rgbScale;

        return float3(r, g, b);
    }

  private:
    bool ownBuffer = false;

  public:
    uint *pixels = 0;
    int width = 0, height = 0;
};
} // namespace Tmpl8