#include "cachesim.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>

uint64_t timestamp = 0;

double sc_time_stamp() { 
  return timestamp;
}

CacheSim::CacheSim() {
  // force random values for uninitialized signals  
  Verilated::randReset(2);

  ram_ = nullptr;
  cache_ = new VVX_cache();

  dram_rsp_active_ = false;
  snp_req_active_ = false;

//#ifdef VCD_OUTPUT
  Verilated::traceEverOn(true);
  trace_ = new VerilatedVcdC;
  cache_->trace(trace_, 99);
  trace_->open("trace.vcd");
//#endif
}

CacheSim::~CacheSim() {
//#ifdef VCD_OUTPUT
  trace_->close();
//#endif
  delete cache_;
  //need to delete the req and rsp vectors
}

void CacheSim::attach_ram(RAM* ram) {
  ram_ = ram;
  dram_rsp_vec_.clear();
}

void CacheSim::reset() {
#ifndef NDEBUG
  std::cout << timestamp << ": [sim] reset()" << std::endl;
#endif

  cache_->reset = 1;
  this->step();
  cache_->reset = 0;
  this->step();

  dram_rsp_vec_.clear();
  //clear req and rsp vecs
  
}

void CacheSim::step() {
  //toggle clock
  cache_->clk = 0;
  this->eval();

  cache_->clk = 1;
  this->eval();

  //handle core and dram reqs and rsps
  this->eval_reqs();
  this->eval_rsps();
  this->eval_dram_bus();
}

void CacheSim::eval() {
  cache_->eval();
//#ifdef VCD_OUTPUT
  trace_->dump(timestamp);
//#endif
  ++timestamp;
}

void CacheSim::run(){
#ifndef NDEBUG
  std::cout << timestamp << ": [sim] run()" << std::endl;
#endif
  this->step();

  int valid = 300; 
  
  while (valid > -1) {
      this->step();
      if(!cache_->core_req_valid && !cache_->core_rsp_valid){
        valid--; 
      }
  }
}

void CacheSim::clear_req(){
  cache_->core_req_valid = 0; 
}


void CacheSim::send_req(core_req_t *req){
  core_req_vec_.push(req);
  unsigned int *data = new unsigned int[4];
  core_rsp_vec_.insert(std::pair<unsigned int, unsigned int*>(req->tag, data));
}

bool CacheSim::get_core_req_ready(){
  return cache_->core_req_ready; 
}

bool CacheSim::get_core_rsp_ready(){
  return cache_->core_rsp_ready; 
}

void CacheSim::eval_reqs(){
  //check to see if cache is accepting reqs
  if(!core_req_vec_.empty() && cache_->core_req_ready){
    core_req_t *req = core_req_vec_.front();

    cache_->core_req_valid = req->valid;
    cache_->core_req_rw = req->rw; 
    cache_->core_req_byteen = req->byteen;

    cache_->core_req_addr[0] = req->addr[0];
    cache_->core_req_addr[1] = req->addr[1];
    cache_->core_req_addr[2] = req->addr[2];
    cache_->core_req_addr[3] = req->addr[3];

    cache_->core_req_data[0] = req->data[0];
    cache_->core_req_data[1] = req->data[1];
    cache_->core_req_data[2] = req->data[2];
    cache_->core_req_data[3] = req->data[3];

    cache_->core_req_tag = req->tag; 

    core_req_vec_.pop();
   
  } else {
    clear_req();
  }
}

void CacheSim::eval_rsps(){
  //check to see if a request has been responded to
  if (cache_->core_rsp_valid){
      core_rsp_vec_.at(cache_->core_rsp_tag)[0] = cache_->core_rsp_data[0];
      core_rsp_vec_.at(cache_->core_rsp_tag)[1] = cache_->core_rsp_data[1];
      core_rsp_vec_.at(cache_->core_rsp_tag)[2] = cache_->core_rsp_data[2];
      core_rsp_vec_.at(cache_->core_rsp_tag)[3] = cache_->core_rsp_data[3];
  }
}

void CacheSim::eval_dram_bus() {
  if (ram_ == nullptr) {
    cache_->dram_req_ready = 0;
    return;
  }

  // schedule DRAM responses
  int dequeue_index = -1;
  for (int i = 0; i < dram_rsp_vec_.size(); i++) {
    if (dram_rsp_vec_[i].cycles_left > 0) {
      dram_rsp_vec_[i].cycles_left -= 1;
    }
    if ((dequeue_index == -1) 
     && (dram_rsp_vec_[i].cycles_left == 0)) {
      dequeue_index = i;
    }
  }

  // send DRAM response  
  if (dram_rsp_active_
   && cache_->dram_rsp_valid 
   && cache_->dram_rsp_ready) {
    dram_rsp_active_ = false;
  }
  if (!dram_rsp_active_) {
    if (dequeue_index != -1) { //time to respond to the request
      cache_->dram_rsp_valid = 1;

      //copy data from the rsp queue to the cache module
      memcpy((uint8_t*)cache_->dram_rsp_data, dram_rsp_vec_[dequeue_index].data, GLOBAL_BLOCK_SIZE);

      cache_->dram_rsp_tag = dram_rsp_vec_[dequeue_index].tag;    
      free(dram_rsp_vec_[dequeue_index].data); //take data out of the queue
      dram_rsp_vec_.erase(dram_rsp_vec_.begin() + dequeue_index);
      dram_rsp_active_ = true;
    } else {
      cache_->dram_rsp_valid = 0;
    }
  }

  // handle DRAM stalls
  bool dram_stalled = false;
#ifdef ENABLE_DRAM_STALLS
  if (0 == ((timestamp/2) % DRAM_STALLS_MODULO)) { 
    dram_stalled = true;
  } else
  if (dram_rsp_vec_.size() >= DRAM_RQ_SIZE) {
    dram_stalled = true;
  }
#endif

  // process DRAM requests
  if (!dram_stalled) {
    if (cache_->dram_req_valid) {
      if (cache_->dram_req_rw) { //write = 1
        uint64_t byteen = cache_->dram_req_byteen;
        unsigned base_addr = (cache_->dram_req_addr * GLOBAL_BLOCK_SIZE);
        uint8_t* data = (uint8_t*)(cache_->dram_req_data);
        for (int i = 0; i < GLOBAL_BLOCK_SIZE; i++) {
          if ((byteen >> i) & 0x1) {            
            (*ram_)[base_addr + i] = data[i];
          }
        }
      } else {
        dram_req_t dram_req;
        dram_req.cycles_left = DRAM_LATENCY;     
        dram_req.data = (uint8_t*)malloc(GLOBAL_BLOCK_SIZE);
        dram_req.tag = cache_->dram_req_tag;
        ram_->read(cache_->dram_req_addr * GLOBAL_BLOCK_SIZE, GLOBAL_BLOCK_SIZE, dram_req.data);
        dram_rsp_vec_.push_back(dram_req);
      } 
    }    
  }

  cache_->dram_req_ready = ~dram_stalled;
}

bool CacheSim::assert_equal(unsigned int* data, unsigned int tag){
  int check = 0;
  unsigned int *rsp = core_rsp_vec_.at(tag);
  for (int i = 0; i < 4; ++i){
    for (int j = 0; j < 4; ++j){
      if (data[i] == rsp[j]){
        check++;
      }
    }
  }

  return check; 

}

//DEBUG

void CacheSim::get_core_rsp(unsigned int (&rsp)[4]){
  rsp[0] = cache_->core_rsp_data[0];
  rsp[1] = cache_->core_rsp_data[1];
  rsp[2] = cache_->core_rsp_data[2];
  rsp[3] = cache_->core_rsp_data[3];
  //std::cout << std::hex << "core_rsp_valid: " << cache_->core_rsp_valid << std::endl;
  //std::cout << std::hex << "core_rsp_data: " << cache_->core_rsp_data << std::endl;
  //std::cout << std::hex << "core_rsp_tag: " << cache_->core_rsp_tag << std::endl; 
}

void CacheSim::get_core_req(){
  char check = cache_->core_req_valid;
  std::cout << std::hex << "core_req_valid: " << check << std::endl;
  std::cout << std::hex << "core_req_data[0]: " << cache_->core_req_data[0] << std::endl;
  std::cout << std::hex << "core_req_data[1]: " << cache_->core_req_data[1] << std::endl;
  std::cout << std::hex << "core_req_data[2]: " << cache_->core_req_data[2] << std::endl;
  std::cout << std::hex << "core_req_data[3]: " << cache_->core_req_data[3] << std::endl;
  std::cout << std::hex << "core_req_tag: " << cache_->core_req_tag << std::endl; 
}

void CacheSim::get_dram_req(){
  std::cout << std::hex << "dram_req_valid: " << cache_->dram_req_valid << std::endl;
  std::cout << std::hex << "dram_req_rw: " << cache_->dram_req_rw << std::endl;
  std::cout << std::hex << "dram_req_byteen: " << cache_->dram_req_byteen << std::endl;
  std::cout << std::hex << "dram_req_addr: " << cache_->dram_req_addr << std::endl;
  std::cout << std::hex << "dram_req_data: " << cache_->dram_req_data << std::endl; 
  std::cout << std::hex << "dram_req_tag: " << cache_->dram_req_tag << std::endl;
}

void CacheSim::get_dram_rsp(){
  std::cout << std::hex << "dram_rsp_valid: " << cache_->dram_rsp_valid << std::endl;
  std::cout << std::hex << "dram_rsp_data: " << cache_->dram_rsp_data << std::endl; 
  std::cout << std::hex << "dram_rsp_tag: " << cache_->dram_rsp_tag << std::endl;
  std::cout << std::hex << "dram_rsp_ready: " << cache_->dram_rsp_ready << std::endl;
}

