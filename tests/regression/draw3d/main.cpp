#include <iostream>
#include <vector>
#include <unistd.h>
#include <string.h>
#include <chrono>
#include <cmath>
#include <array>
#include <assert.h>
#include <vortex.h>
#include "common.h"
#include "utils.h"
#include <cocogfx/include/cgltrace.hpp>

using namespace cocogfx;

#define RT_CHECK(_expr)                                         \
   do {                                                         \
     int _ret = _expr;                                          \
     if (0 == _ret)                                             \
       break;                                                   \
     printf("Error: '%s' returned %d!\n", #_expr, (int)_ret);   \
	 cleanup();			                                              \
     exit(-1);                                                  \
   } while (false)

///////////////////////////////////////////////////////////////////////////////

const char* kernel_file = "kernel.bin";
const char* trace_file  = "triangle.cgltrace";
const char* output_file = "output.png";
const char* reference_file  = nullptr;

uint32_t clear_color = 0xFFFFFFFF;
uint32_t clear_depth = 0xFFFFFFFF;

uint32_t dst_width  = 128;
uint32_t dst_height = 128;

uint32_t zbuf_stride = 4;
uint32_t zbuf_pitch  = dst_width * zbuf_stride;
uint32_t zbuf_size   = dst_height * zbuf_pitch;

uint32_t cbuf_stride = 4;
uint32_t cbuf_pitch  = dst_width * cbuf_stride;
uint32_t cbuf_size   = dst_width * cbuf_pitch;

vx_device_h device = nullptr;
vx_buffer_h staging_buf = nullptr;

uint64_t zbuf_addr = -1;
uint64_t cbuf_addr = -1;
uint64_t texbuf_addr = -1;
uint64_t tilebuf_addr = -1;
uint64_t primbuf_addr = -1;

kernel_arg_t kernel_arg;

uint32_t tile_size = 1 << RASTER_TILE_LOGSIZE;

static void show_usage() {
   std::cout << "Vortex 3D Rendering Test." << std::endl;
   std::cout << "Usage: [-t trace] [-o output] [-r reference] [-w width] [-h height]" << std::endl;
}

static void parse_args(int argc, char **argv) {
  int c;
  while ((c = getopt(argc, argv, "t:i:o:r:w:h:t:?")) != -1) {
    switch (c) {
    case 't':
      trace_file = optarg;
      break;
    case 'o':
      output_file = optarg;
      break;
    case 'r':
      reference_file = optarg;
      break;
    case 'w':
      dst_width = std::atoi(optarg);
      break;
    case 'h':
      dst_height = std::atoi(optarg);
      break;
    case '?': {
      show_usage();
      exit(0);
    } break;
    default:
      show_usage();
      exit(-1);
    }
  }
}

void cleanup() {
  if (staging_buf) {
    vx_buf_free(staging_buf);
  }
  if (device) {     
    if (zbuf_addr != -1ull) vx_mem_free(device, zbuf_addr);
    if (cbuf_addr != -1ull) vx_mem_free(device, cbuf_addr);
    if (texbuf_addr != -1ull) vx_mem_free(device, texbuf_addr);
    if (tilebuf_addr != -1ull) vx_mem_free(device, tilebuf_addr);
    if (primbuf_addr != -1ull) vx_mem_free(device, primbuf_addr);
    vx_dev_close(device);
  }
}

int render(const CGLTrace& trace) {
  // render each draw call
  for (auto& drawcall : trace.drawcalls) {
    auto& states = drawcall.states;

    if (states.texture_enabled) {
      std::vector<uint8_t> texbuf;    
      std::vector<uint32_t> mip_offsets;

      auto& texture = trace.textures.at(drawcall.texture_id);      
      
      auto tex_bpp = Format::GetInfo(texture.format).BytePerPixel;
      auto tex_pitch = texture.width * tex_bpp;

      // generate mipmaps
      RT_CHECK(GenerateMipmaps(texbuf, mip_offsets, texture.pixels.data(), texture.format, texture.width, texture.height, tex_pitch));

      uint32_t tex_logwidth = log2ceil(texture.width);
      uint32_t tex_logheight = log2ceil(texture.height);

      int tex_format = toVXFormat(texture.format);

      int tex_filter = (states.texture_magfilter != CGLTrace::FILTER_NEAREST) 
                    || (states.texture_magfilter != CGLTrace::FILTER_NEAREST);

      int tex_wrapU = (states.texture_addressU == CGLTrace::ADDRESS_WRAP);
      int tex_wrapV = (states.texture_addressU == CGLTrace::ADDRESS_WRAP);

      // allocate texture memory
      if (texbuf_addr != -1ull) vx_mem_free(device, texbuf_addr); 
      RT_CHECK(vx_mem_alloc(device, texbuf.size(), &texbuf_addr));
      std::cout << "texbuf_addr=0x" << std::hex << texbuf_addr << std::endl;

      // upload texture data
      std::cout << "upload texture buffer" << std::endl;      
      {    
        RT_CHECK(vx_buf_alloc(device, texbuf.size(), &staging_buf));
        auto buf_ptr = (uint8_t*)vx_host_ptr(staging_buf);
        memcpy(buf_ptr, texbuf.data(), texbuf.size());
        RT_CHECK(vx_copy_to_dev(staging_buf, texbuf_addr, texbuf.size(), 0));
        vx_buf_free(staging_buf);
        staging_buf = nullptr;
      }

      // configure texture units
      vx_dcr_write(device, DCR_TEX_STAGE,  0);
      vx_dcr_write(device, DCR_TEX_LOGDIM, (tex_logheight << 16) | tex_logwidth);	
      vx_dcr_write(device, DCR_TEX_FORMAT, tex_format);
      vx_dcr_write(device, DCR_TEX_WRAP,   (tex_wrapV << 16) | tex_wrapU);
      vx_dcr_write(device, DCR_TEX_FILTER, tex_filter);
      vx_dcr_write(device, DCR_TEX_ADDR,   texbuf_addr);
      for (uint32_t i = 0; i < mip_offsets.size(); ++i) {
        assert(i < TEX_LOD_MAX);
        vx_dcr_write(device, DCR_TEX_MIPOFF(i), mip_offsets.at(i));
      };
    }

    std::vector<uint8_t> tilebuf;
    std::vector<uint8_t> primbuf;
    
    // Perform tile binning
    auto num_tiles = Binning(tilebuf, primbuf, drawcall.vertices, drawcall.primitives, dst_width, dst_height, drawcall.viewport.near, drawcall.viewport.far, tile_size);
    std::cout << "Binning allocated " << num_tiles << " tiles and " << primbuf.size() << " primitives." << std::endl;

    // allocate tile memory
    if (tilebuf_addr != -1ull) vx_mem_free(device, tilebuf_addr); 
    if (primbuf_addr != -1ull) vx_mem_free(device, primbuf_addr); 
    RT_CHECK(vx_mem_alloc(device, tilebuf.size(), &tilebuf_addr));
    RT_CHECK(vx_mem_alloc(device, primbuf.size(), &primbuf_addr));
    std::cout << "tilebuf_addr=0x" << std::hex << tilebuf_addr << std::endl;
    std::cout << "primbuf_addr=0x" << std::hex << primbuf_addr << std::endl;

    uint32_t alloc_size = std::max({tilebuf.size(), primbuf.size(), sizeof(kernel_arg_t)});
    RT_CHECK(vx_buf_alloc(device, alloc_size, &staging_buf));
    
    // upload tiles buffer
    std::cout << "upload tiles buffer" << std::endl;      
    {    
      auto buf_ptr = (uint8_t*)vx_host_ptr(staging_buf);
      memcpy(buf_ptr, tilebuf.data(), tilebuf.size());
      RT_CHECK(vx_copy_to_dev(staging_buf, tilebuf_addr, tilebuf.size(), 0));
    }

    // upload primitives buffer
    std::cout << "upload primitives buffer" << std::endl;      
    {    
      auto buf_ptr = (uint8_t*)vx_host_ptr(staging_buf);
      memcpy(buf_ptr, primbuf.data(), primbuf.size());
      RT_CHECK(vx_copy_to_dev(staging_buf, primbuf_addr, primbuf.size(), 0));
    }

    // upload kernel argument
    std::cout << "upload kernel argument" << std::endl;
    {
      kernel_arg.depth_enabled = states.depth_test;
      kernel_arg.color_enabled = states.color_enabled;
      kernel_arg.tex_enabled   = states.texture_enabled;
      kernel_arg.tex_modulate  = (states.texture_enabled && states.texture_envmode == CGLTrace::ENVMODE_MODULATE);
      kernel_arg.prim_addr     = primbuf_addr;

      if (kernel_arg.tex_modulate)
        kernel_arg.color_enabled = false;
      
      auto buf_ptr = (uint8_t*)vx_host_ptr(staging_buf);
      memcpy(buf_ptr, &kernel_arg, sizeof(kernel_arg_t));
      RT_CHECK(vx_copy_to_dev(staging_buf, KERNEL_ARG_DEV_MEM_ADDR, sizeof(kernel_arg_t), 0));
    }

    vx_buf_free(staging_buf);
    staging_buf = nullptr;

    uint32_t primbuf_stride = sizeof(rast_prim_t);

    // configure raster units
    vx_dcr_write(device, DCR_RASTER_TBUF_ADDR,   tilebuf_addr);
    vx_dcr_write(device, DCR_RASTER_TILE_COUNT,  num_tiles);
    vx_dcr_write(device, DCR_RASTER_PBUF_ADDR,   primbuf_addr);
    vx_dcr_write(device, DCR_RASTER_PBUF_STRIDE, primbuf_stride);

    // configure rop color buffer
    vx_dcr_write(device, DCR_ROP_CBUF_ADDR,  cbuf_addr);
    vx_dcr_write(device, DCR_ROP_CBUF_PITCH, cbuf_pitch);
    vx_dcr_write(device, DCR_ROP_CBUF_MASK, states.color_writemask);

    if (states.depth_test || states.stencil_test) {
      // configure rop depth buffer
      vx_dcr_write(device, DCR_ROP_ZBUF_ADDR,  zbuf_addr);
      vx_dcr_write(device, DCR_ROP_ZBUF_PITCH, zbuf_pitch);    
    }

    if (states.depth_test) {    
      // configure rop depth states
      auto depth_func = toVXCompare(states.depth_func);
      vx_dcr_write(device, DCR_ROP_DEPTH_FUNC, depth_func);
      vx_dcr_write(device, DCR_ROP_DEPTH_MASK, states.depth_writemask);
    }

    if (states.stencil_test) {
      // configure rop stencil states
      auto stencil_func  = toVXCompare(states.stencil_func);
      auto stencil_zpass = toVXStencilOp(states.stencil_zpass);
      auto stencil_zfail = toVXStencilOp(states.stencil_zfail);
      auto stencil_fail  = toVXStencilOp(states.stencil_fail);
      vx_dcr_write(device, DCR_ROP_STENCIL_FUNC, stencil_func);
      vx_dcr_write(device, DCR_ROP_STENCIL_ZPASS, stencil_zpass);
      vx_dcr_write(device, DCR_ROP_STENCIL_ZPASS, stencil_zfail);
      vx_dcr_write(device, DCR_ROP_STENCIL_FAIL, stencil_fail);
      vx_dcr_write(device, DCR_ROP_STENCIL_MASK, states.stencil_mask);
      vx_dcr_write(device, DCR_ROP_STENCIL_REF, states.stencil_ref);
    }

    if (states.blend_enabled) {
      // configure rop blend states
      auto blend_src = toVXBlendFunc(states.blend_src);
      auto blend_dst = toVXBlendFunc(states.blend_dst);
      vx_dcr_write(device, DCR_ROP_BLEND_MODE, (ROP_BLEND_MODE_ADD << 16)  // DST
                                             | (ROP_BLEND_MODE_ADD << 0)); // SRC
      vx_dcr_write(device, DCR_ROP_BLEND_FUNC, (blend_dst << 24)  // DST_A
                                             | (blend_dst << 16)  // DST_RGB 
                                             | (blend_src << 8)   // SRC_A
                                             | (blend_src << 0)); // SRC_RGB
    }

    auto time_start = std::chrono::high_resolution_clock::now();

    // start device
    std::cout << "start device" << std::endl;
    RT_CHECK(vx_start(device));

    // wait for completion
    std::cout << "wait for completion" << std::endl;
    RT_CHECK(vx_ready_wait(device, MAX_TIMEOUT));
    
    auto time_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    printf("Elapsed time: %lg ms\n", elapsed);

  }

  // download destination buffer
  std::vector<uint8_t> dst_pixels(cbuf_size);
  {
    std::cout << "download destination buffer" << std::endl;
    RT_CHECK(vx_buf_alloc(device, cbuf_size, &staging_buf));
    RT_CHECK(vx_copy_from_dev(staging_buf, cbuf_addr, cbuf_size, 0));    
    auto buf_ptr = (uint8_t*)vx_host_ptr(staging_buf);
    memcpy(dst_pixels.data(), buf_ptr, cbuf_size);
    vx_buf_free(staging_buf);
    staging_buf = nullptr;
  }

  // save output image
  std::cout << "save output image" << std::endl;
  {
    // save image upside down
    auto bits = dst_pixels.data() + (dst_height-1) * cbuf_pitch;
    RT_CHECK(SaveImage(output_file, FORMAT_A8R8G8B8, bits, dst_width, dst_height, -cbuf_pitch));
  }

  return 0;
}

int main(int argc, char *argv[]) {  
  // parse command arguments
  parse_args(argc, argv);

  if (!ispow2(dst_width)) {
    std::cout << "Error: only power of two dst_width supported: dst_width=" << dst_width << std::endl;    
    return -1;
  }

  if (!ispow2(dst_height)) {
    std::cout << "Error: only power of two dst_height supported: dst_height=" << dst_height << std::endl;
    return -1;
  }

  if (0 != (dst_width % tile_size)) {
    std::cout << "Error: dst_with must be divisible by tile_size" << std::endl;
    return -1;
  }

  if (0 != (dst_height % tile_size)) {
    std::cout << "Error: dst_height must be divisible by tile_size" << std::endl;
    return -1;
  }

  // open device connection
  std::cout << "open device connection" << std::endl;  
  RT_CHECK(vx_dev_open(&device));

  uint64_t isa_flags;
  RT_CHECK(vx_dev_caps(device, VX_CAPS_ISA_FLAGS, &isa_flags));
  if (0 == (isa_flags & (VX_ISA_EXT_RASTER | VX_ISA_EXT_ROP))) {
    std::cout << "raster or rop extensions not supported!" << std::endl;
    cleanup();
    return -1;
  }

  CGLTrace trace;    
  RT_CHECK(trace.load(trace_file));

  // upload program
  std::cout << "upload program" << std::endl;  
  RT_CHECK(vx_upload_kernel_file(device, kernel_file));

  // allocate device memory  
  RT_CHECK(vx_mem_alloc(device, zbuf_size, &zbuf_addr));
  RT_CHECK(vx_mem_alloc(device, cbuf_size, &cbuf_addr));

  std::cout << "zbuf_addr=0x" << std::hex << zbuf_addr << std::endl;
  std::cout << "cbuf_addr=0x" << std::hex << cbuf_addr << std::endl;

  // allocate staging buffer  
  std::cout << "allocate staging buffer" << std::endl;    
  uint32_t alloc_size = std::max(zbuf_size, cbuf_size);
  RT_CHECK(vx_buf_alloc(device, alloc_size, &staging_buf));
  
  // clear depth buffer
  std::cout << "clear depth buffer" << std::endl;      
  {    
    auto buf_ptr = (uint32_t*)vx_host_ptr(staging_buf);
    for (uint32_t i = 0; i < (zbuf_size/4); ++i) {
      buf_ptr[i] = clear_depth;
    }    
    RT_CHECK(vx_copy_to_dev(staging_buf, zbuf_addr, zbuf_size, 0));  
  }

  // clear destination buffer
  std::cout << "clear destination buffer" << std::endl;      
  {    
    auto buf_ptr = (uint32_t*)vx_host_ptr(staging_buf);
    for (uint32_t i = 0; i < (cbuf_size/4); ++i) {
      buf_ptr[i] = clear_color;
    }    
    RT_CHECK(vx_copy_to_dev(staging_buf, cbuf_addr, cbuf_size, 0));  
  }  

  vx_buf_free(staging_buf);
  staging_buf = nullptr;

  // update kernel arguments
  kernel_arg.dst_width     = dst_width;
  kernel_arg.dst_height    = dst_height;

  kernel_arg.cbuf_stride   = cbuf_stride;
  kernel_arg.cbuf_pitch    = cbuf_pitch;    
  kernel_arg.cbuf_addr     = cbuf_addr;

  kernel_arg.zbuf_stride   = zbuf_stride;
  kernel_arg.zbuf_pitch    = zbuf_pitch;    
  kernel_arg.zbuf_addr     = zbuf_addr;

  // run tests
  std::cout << "render" << std::endl;
  RT_CHECK(render(trace));

  // cleanup
  std::cout << "cleanup" << std::endl;  
  cleanup();  

  if (reference_file) {
    auto errors = CompareImages(output_file, reference_file, FORMAT_A8R8G8B8);
    if (0 == errors) {
      std::cout << "PASSED!" << std::endl;
    } else {
      std::cout << "FAILED!" << std::endl;
      return errors;
    }
  } 

  return 0;
}