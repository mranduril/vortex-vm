#include <cstdint>
#include <vector>
#include <bitmanip.h>
#include <cocogfx/include/format.hpp>
#include <cocogfx/include/blitter.hpp>

int LoadImage(const char *filename,
              cocogfx::ePixelFormat format,
              std::vector<uint8_t> &pixels,
              uint32_t *width,
              uint32_t *height);

int SaveImage(const char *filename,
              cocogfx::ePixelFormat format,
              const uint8_t* pixels,
              uint32_t width,
              uint32_t height,
              int32_t pitch);

void dump_image(const std::vector<uint8_t>& pixels, 
                uint32_t width, 
                uint32_t height, 
                uint32_t bpp);

int CompareImages(const char* filename1, 
                  const char* filename2, 
                  cocogfx::ePixelFormat format);