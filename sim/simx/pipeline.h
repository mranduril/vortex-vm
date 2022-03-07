
#pragma once

#include <memory>
#include <iostream>
#include <util.h>
#include "types.h"
#include "arch.h"
#include "debug.h"

namespace vortex {

class ITraceData {
public:
    ITraceData() {}
    virtual ~ITraceData() {}
};

struct LsuTraceData : public ITraceData {
  std::vector<mem_addr_size_t> mem_addrs;
  LsuTraceData(uint32_t num_threads) : mem_addrs(num_threads) {}
};

struct GPUTraceData : public ITraceData {
  const WarpMask active_warps;
  GPUTraceData(const WarpMask& active_warps) : active_warps(active_warps) {}
};

struct pipeline_trace_t {
public:
  //--
  const uint64_t  uuid;
  
  //--
  uint32_t    cid;
  uint32_t    wid;  
  ThreadMask  tmask;
  Word        PC;

  //--
  uint32_t    rdest;
  RegType     rdest_type;
  bool        wb;

  //--
  RegMask     used_iregs;
  RegMask     used_fregs;
  RegMask     used_vregs;

  //- 
  ExeType     exe_type; 

  //--
  union {
    LsuType lsu_type;
    AluType alu_type;
    FpuType fpu_type;
    GpuType gpu_type;
  };

  ITraceData* data;

  bool fetch_stall;

private:
  bool stalled_;

public:
  pipeline_trace_t(uint64_t uuid) 
    : uuid(uuid)
    , cid(0)
    , wid(0)
    , PC(0)    
    , rdest(0)
    , rdest_type(RegType::None)
    , wb(false)
    , exe_type(ExeType::NOP)
    , data(nullptr)
    , fetch_stall(false)
    , stalled_(false) 
  {}
  
  ~pipeline_trace_t() {
    if (data)
      delete data;
  }

  bool suspend() {
    bool old = stalled_;
    stalled_ = true;
    return old;
  }

  void resume() {
    stalled_ = false;
  }
};

inline std::ostream &operator<<(std::ostream &os, const pipeline_trace_t& state) {
  os << "coreid=" << state.cid << ", wid=" << state.wid << ", PC=" << std::hex << state.PC;
  os << ", wb=" << state.wb;
  if (state.wb) {
     os << ", rd=" << state.rdest_type << std::dec << state.rdest;
  }
  os << ", ex=" << state.exe_type;
  os << " (#" << std::dec << state.uuid << ")";
  return os;
}

class PipelineLatch {
protected:
  const char* name_;
  std::queue<pipeline_trace_t*> queue_;

public:
  PipelineLatch(const char* name = nullptr) 
    : name_(name) 
  {}
  
  bool empty() const {
    return queue_.empty();
  }

  pipeline_trace_t* front() {
    return queue_.front();
  }

  pipeline_trace_t* back() {
    return queue_.back();
  }

  void push(pipeline_trace_t* value) {    
    queue_.push(value);
  }

  void pop() {
    queue_.pop();
  }

  void clear() {
    std::queue<pipeline_trace_t*> empty;
    std::swap(queue_, empty );
  }
};

}