#pragma once

#include <stdint.h>

namespace vortex {

class Arch;
class RAM;
class ProcessorImpl;

class Processor {
public:
  Processor(const Arch& arch);
  ~Processor();

  void attach_ram(RAM* mem);

  int run();

  void write_dcr(uint32_t addr, uint32_t value);

  uint32_t get_satp();
  void set_satp(uint32_t satp);

private:
  ProcessorImpl* impl_;
  uint32_t satp;
};

}
