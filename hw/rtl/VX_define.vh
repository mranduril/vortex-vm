`ifndef VX_DEFINE
`define VX_DEFINE

`include "VX_platform.vh"
`include "VX_config.vh"
`include "VX_types.vh"

///////////////////////////////////////////////////////////////////////////////

`define NW_BITS         `LOG2UP(`NUM_WARPS)

`define NT_BITS         `LOG2UP(`NUM_THREADS)

`define NC_BITS         `LOG2UP(`NUM_CORES)

`define NB_BITS         `LOG2UP(`NUM_BARRIERS)

`define NUM_IREGS       32

`define NRI_BITS        `LOG2UP(`NUM_IREGS)

`ifdef EXT_F_ENABLE
`define NUM_REGS        (2 * `NUM_IREGS)
`else
`define NUM_REGS        `NUM_IREGS
`endif

`define NR_BITS         `LOG2UP(`NUM_REGS)

`define CSR_ADDR_BITS   12

`define DCR_ADDR_BITS   12

`define PERF_CTR_BITS   44

`define UUID_BITS       44

///////////////////////////////////////////////////////////////////////////////

`define EX_NOP          3'h0
`define EX_ALU          3'h1
`define EX_LSU          3'h2
`define EX_CSR          3'h3
`define EX_FPU          3'h4
`define EX_GPU          3'h5
`define EX_BITS         3

///////////////////////////////////////////////////////////////////////////////

`define INST_LUI        7'b0110111
`define INST_AUIPC      7'b0010111
`define INST_JAL        7'b1101111
`define INST_JALR       7'b1100111
`define INST_B          7'b1100011 // branch instructions
`define INST_L          7'b0000011 // load instructions
`define INST_S          7'b0100011 // store instructions
`define INST_I          7'b0010011 // immediate instructions
`define INST_R          7'b0110011 // register instructions
`define INST_FENCE      7'b0001111 // Fence instructions
`define INST_SYS        7'b1110011 // system instructions

`define INST_FL         7'b0000111 // float load instruction
`define INST_FS         7'b0100111 // float store  instruction
`define INST_FMADD      7'b1000011  
`define INST_FMSUB      7'b1000111
`define INST_FNMSUB     7'b1001011
`define INST_FNMADD     7'b1001111 
`define INST_FCI        7'b1010011 // float common instructions

// Custom extension opcodes
`define INST_EXT1       7'b0001011 // 0x0B
`define INST_EXT2       7'b0101011 // 0x2B
`define INST_EXT3       7'b1011011 // 0x5B
`define INST_EXT4       7'b1111011 // 0x7B

///////////////////////////////////////////////////////////////////////////////

`define INST_FRM_RNE    3'b000  // round to nearest even
`define INST_FRM_RTZ    3'b001  // round to zero
`define INST_FRM_RDN    3'b010  // round to -inf
`define INST_FRM_RUP    3'b011  // round to +inf
`define INST_FRM_RMM    3'b100  // round to nearest max magnitude
`define INST_FRM_DYN    3'b111  // dynamic mode
`define INST_FRM_BITS   3

///////////////////////////////////////////////////////////////////////////////

`define INST_OP_BITS    4
`define INST_MOD_BITS   3

///////////////////////////////////////////////////////////////////////////////

`define INST_ALU_ADD         4'b0000
`define INST_ALU_LUI         4'b0010
`define INST_ALU_AUIPC       4'b0011
`define INST_ALU_SLTU        4'b0100
`define INST_ALU_SLT         4'b0101
`define INST_ALU_SRL         4'b1000
`define INST_ALU_SRA         4'b1001
`define INST_ALU_SUB         4'b1011
`define INST_ALU_AND         4'b1100
`define INST_ALU_OR          4'b1101
`define INST_ALU_XOR         4'b1110
`define INST_ALU_SLL         4'b1111
`define INST_ALU_OTHER       4'b0111
`define INST_ALU_BITS        4
`define INST_ALU_OP(x)       x[`INST_ALU_BITS-1:0]
`define INST_ALU_OP_CLASS(x) x[3:2]
`define INST_ALU_SIGNED(x)   x[0]
`define INST_ALU_IS_BR(x)    x[0]
`define INST_ALU_IS_MUL(x)   x[1]

`define INST_BR_EQ           4'b0000
`define INST_BR_NE           4'b0010
`define INST_BR_LTU          4'b0100 
`define INST_BR_GEU          4'b0110 
`define INST_BR_LT           4'b0101
`define INST_BR_GE           4'b0111
`define INST_BR_JAL          4'b1000
`define INST_BR_JALR         4'b1001
`define INST_BR_ECALL        4'b1010
`define INST_BR_EBREAK       4'b1011
`define INST_BR_URET         4'b1100
`define INST_BR_SRET         4'b1101
`define INST_BR_MRET         4'b1110
`define INST_BR_OTHER        4'b1111
`define INST_BR_BITS         4
`define INST_BR_NEG(x)       x[1]
`define INST_BR_LESS(x)      x[2]
`define INST_BR_STATIC(x)    x[3]

`define INST_MUL_MUL         3'h0
`define INST_MUL_MULH        3'h1
`define INST_MUL_MULHSU      3'h2
`define INST_MUL_MULHU       3'h3
`define INST_MUL_DIV         3'h4
`define INST_MUL_DIVU        3'h5
`define INST_MUL_REM         3'h6
`define INST_MUL_REMU        3'h7
`define INST_MUL_BITS        3
`define INST_MUL_IS_DIV(x)   x[2]

`define INST_FMT_B           3'b000
`define INST_FMT_H           3'b001
`define INST_FMT_W           3'b010
`define INST_FMT_BU          3'b100
`define INST_FMT_HU          3'b101

`define INST_LSU_LB          4'b0000 
`define INST_LSU_LH          4'b0001
`define INST_LSU_LW          4'b0010
`define INST_LSU_LBU         4'b0100
`define INST_LSU_LHU         4'b0101
`define INST_LSU_SB          4'b1000 
`define INST_LSU_SH          4'b1001
`define INST_LSU_SW          4'b1010
`define INST_LSU_BITS        4
`define INST_LSU_FMT(x)      x[2:0]
`define INST_LSU_WSIZE(x)    x[1:0]
`define INST_LSU_IS_MEM(x)   (3'h0 == x)
`define INST_LSU_IS_FENCE(x) (3'h1 == x)
`define INST_LSU_IS_PREFETCH(x) (3'h2 == x)

`define INST_FENCE_BITS      1
`define INST_FENCE_D         1'h0
`define INST_FENCE_I         1'h1

`define INST_CSR_RW          2'h1
`define INST_CSR_RS          2'h2
`define INST_CSR_RC          2'h3
`define INST_CSR_OTHER       2'h0
`define INST_CSR_BITS        2

`define INST_FPU_ADD         4'h0 
`define INST_FPU_SUB         4'h4 
`define INST_FPU_MUL         4'h8 
`define INST_FPU_DIV         4'hC
`define INST_FPU_CVTWS       4'h1  // FCVT.W.S
`define INST_FPU_CVTWUS      4'h5  // FCVT.WU.S
`define INST_FPU_CVTSW       4'h9  // FCVT.S.W
`define INST_FPU_CVTSWU      4'hD  // FCVT.S.WU
`define INST_FPU_SQRT        4'h2
`define INST_FPU_CLASS       4'h6  
`define INST_FPU_CMP         4'hA
`define INST_FPU_MISC        4'hE  // SGNJ, SGNJN, SGNJX, FMIN, FMAX, MVXW, MVWX 
`define INST_FPU_MADD        4'h3 
`define INST_FPU_MSUB        4'h7   
`define INST_FPU_NMSUB       4'hB   
`define INST_FPU_NMADD       4'hF
`define INST_FPU_BITS        4

`define INST_GPU_TMC         4'h0
`define INST_GPU_WSPAWN      4'h1 
`define INST_GPU_SPLIT       4'h2
`define INST_GPU_JOIN        4'h3
`define INST_GPU_BAR         4'h4
`define INST_GPU_PRED        4'h5

`define INST_GPU_TEX         4'h6
`define INST_GPU_RASTER      4'h7
`define INST_GPU_ROP         4'h8
`define INST_GPU_CMOV        4'h9
`define INST_GPU_IMADD       4'hA
`define INST_GPU_BITS        4

///////////////////////////////////////////////////////////////////////////////

// non-cacheable tag bits
`define NC_TAG_BIT              1

// texture tag bits
`define TEX_TAG_BIT             1

// cache address type bits
`define CACHE_ADDR_TYPE_BITS    (`NC_TAG_BIT + `SM_ENABLED)

////////////////////////// Icache Configurable Knobs //////////////////////////

// Cache ID
`define ICACHE_ID               $sformatf("core%0d-icache", CORE_ID)

// Word size in bytes
`define ICACHE_WORD_SIZE        4
`define ICACHE_ADDR_WIDTH       (32-`CLOG2(`ICACHE_WORD_SIZE))

// Block size in bytes
`define ICACHE_LINE_SIZE        `L1_BLOCK_SIZE

// Response tag select bits       
`define ICACHE_TAG_SEL_BITS     `NW_BITS

// Core request tag bits
`define ICACHE_TAG_WIDTH        (`UUID_BITS + `ICACHE_TAG_SEL_BITS)

// Input request size
`define ICACHE_NUM_REQS         1

// Memory request data bits
`define ICACHE_MEM_DATA_WIDTH   (`ICACHE_LINE_SIZE * 8)

// Memory request address bits
`define ICACHE_MEM_ADDR_WIDTH   (32 - `CLOG2(`ICACHE_LINE_SIZE))

// Memory request tag bits
`define ICACHE_MEM_TAG_WIDTH    `CLOG2(`ICACHE_MSHR_SIZE)

////////////////////////// Dcache Configurable Knobs //////////////////////////

// Cache ID
`define DCACHE_ID               $sformatf("core%0d-dcache", CORE_ID)

// Word size in bytes
`define DCACHE_WORD_SIZE        4
`define DCACHE_ADDR_WIDTH       (32-`CLOG2(`DCACHE_WORD_SIZE))

// Block size in bytes
`define DCACHE_LINE_SIZE        `L1_BLOCK_SIZE

// Response tag select bits
`define LSUQ_ADDR_BITS          `LOG2UP(`LSUQ_SIZE)
`define DCACHE_TAG_SEL_BITS     (`LSUQ_ADDR_BITS + `CACHE_ADDR_TYPE_BITS)

// Core request tag bits
`define DCACHE_TAG_WIDTH        (`UUID_BITS + `DCACHE_TAG_SEL_BITS)
 
// Memory request data bits
`define DCACHE_MEM_DATA_WIDTH   (`DCACHE_LINE_SIZE * 8)

// Memory request address bits
`define DCACHE_MEM_ADDR_WIDTH   (32 - `CLOG2(`DCACHE_LINE_SIZE))

// Memory byte enable bits
`define DCACHE_MEM_BYTEEN_WIDTH `DCACHE_LINE_SIZE

// Input request size
`define DCACHE_NUM_REQS         `NUM_THREADS

// Memory request tag bits
`define DCACHE_SMEM_TAG_SEL_BITS (`DCACHE_TAG_SEL_BITS - `SM_ENABLED)
`define DCACHE_SMEM_TAG_WIDTH    (`UUID_BITS + `DCACHE_SMEM_TAG_SEL_BITS)
`ifdef EXT_TEX_ENABLE
`define DCACHE_TEX_TAG_SEL_BITS  `MAX(`DCACHE_SMEM_TAG_SEL_BITS, `TCACHE_TAG_SEL_BITS)
`define DCACHE_TEX_TAG_WIDTH     (`UUID_BITS + `DCACHE_TEX_TAG_SEL_BITS)
`else 
`define DCACHE_TEX_TAG_SEL_BITS  `DCACHE_SMEM_TAG_SEL_BITS
`define DCACHE_TEX_TAG_WIDTH     `DCACHE_SMEM_TAG_WIDTH
`endif
`define _DMEM_ADDR_RATIO_W      `CLOG2(`DCACHE_LINE_SIZE / `DCACHE_WORD_SIZE)
`define _DNC_MEM_TAG_WIDTH      (`CLOG2(`DCACHE_NUM_REQS) + `_DMEM_ADDR_RATIO_W + (`DCACHE_TEX_TAG_WIDTH + `EXT_TEX_ENABLED))
`define DCACHE_MEM_TAG_WIDTH    `MAX((`CLOG2(`DCACHE_NUM_BANKS) + `CLOG2(`DCACHE_MSHR_SIZE) + `NC_TAG_BIT), `_DNC_MEM_TAG_WIDTH)

// Merged D-cache/I-cache memory tag
`define L1_MEM_RGB_TAG_WIDTH    `MAX(`ICACHE_MEM_TAG_WIDTH, `DCACHE_MEM_TAG_WIDTH)
`define L1_MEM_TAG_WIDTH        (`L1_MEM_RGB_TAG_WIDTH + `CLOG2(2))

////////////////////////// SM Configurable Knobs //////////////////////////////

// Cache ID
`define SMEM_ID                 $sformatf("core%0d-smem", CORE_ID)

// bank address offset
`define SMEM_BANK_ADDR_OFFSET   `CLOG2(`STACK_SIZE / `DCACHE_WORD_SIZE)

////////////////////////// L2cache Configurable Knobs /////////////////////////

// Cache ID
`define L2_CACHE_ID             $sformatf("cluster%0d-l2cache", CLUSTER_ID)

// Word size in bytes
`define L2_WORD_SIZE            `DCACHE_LINE_SIZE
`define L2_ADDR_WIDTH           (32-`CLOG2(`L2_WORD_SIZE))

// Block size in bytes
`define L2_CACHE_LINE_SIZE      (`L2_ENABLED ? `MEM_BLOCK_SIZE : `L2_WORD_SIZE)

// Memory request data bits
`define L2_MEM_DATA_WIDTH       (`L2_CACHE_LINE_SIZE * 8)

// Memory request address bits
`define L2_MEM_ADDR_WIDTH       (32 - `CLOG2(`L2_CACHE_LINE_SIZE))

// Memory byte enable bits
`define L2_MEM_BYTEEN_WIDTH     `L2_CACHE_LINE_SIZE

// Input request size
`define L2_NUM_REQS             `NUM_CORES

// Memory request tag bits
`define _L2_MEM_ADDR_RATIO_W    `CLOG2(`L2_CACHE_LINE_SIZE / `L2_WORD_SIZE)
`define _L2_NC_MEM_TAG_WIDTH    (`CLOG2(`L2_NUM_REQS) + `_L2_MEM_ADDR_RATIO_W + `L1_MEM_TAG_WIDTH)
`define _L2_MEM_TAG_WIDTH       `MAX((`CLOG2(`L2_NUM_BANKS) + `CLOG2(`L2_MSHR_SIZE) + `NC_TAG_BIT), `_L2_NC_MEM_TAG_WIDTH)
`define L2X_MEM_TAG_WIDTH       (`L2_ENABLED ? `_L2_MEM_TAG_WIDTH : (`L1_MEM_TAG_WIDTH + `CLOG2(`L2_NUM_REQS)))
`define L2_MEM_TAG_WIDTH        (`L2X_MEM_TAG_WIDTH + `CLOG2(1 + `EXT_RASTER_ENABLED + `EXT_ROP_ENABLED))

////////////////////////// L3cache Configurable Knobs /////////////////////////

// Cache ID
`define L3_CACHE_ID             "l3cache"

// Word size in bytes
`define L3_WORD_SIZE            `L2_CACHE_LINE_SIZE
`define L3_ADDR_WIDTH           (32-`CLOG2(`L3_WORD_SIZE))

// Block size in bytes
`define L3_CACHE_LINE_SIZE      (`L3_ENABLED ? `MEM_BLOCK_SIZE : `L3_WORD_SIZE)

// Memory request data bits
`define L3_MEM_DATA_WIDTH       (`L3_CACHE_LINE_SIZE * 8)

// Memory request address bits
`define L3_MEM_ADDR_WIDTH       (32 - `CLOG2(`L3_CACHE_LINE_SIZE))

// Memory byte enable bits
`define L3_MEM_BYTEEN_WIDTH     `L3_CACHE_LINE_SIZE

// Input request size
`define L3_NUM_REQS             `NUM_CLUSTERS

// Memory request tag bits
`define _L3_MEM_ADDR_RATIO_W    `CLOG2(`L3_CACHE_LINE_SIZE / `L3_WORD_SIZE)
`define _L3_NC_MEM_TAG_WIDTH    (`CLOG2(`L3_NUM_REQS) + `_L3_MEM_ADDR_RATIO_W + `L2_MEM_TAG_WIDTH)
`define _L3_MEM_TAG_WIDTH       `MAX((`CLOG2(`L3_NUM_BANKS) + `CLOG2(`L3_MSHR_SIZE) + `NC_TAG_BIT), `_L3_NC_MEM_TAG_WIDTH)
`define L3_MEM_TAG_WIDTH        (`L3_ENABLED ? `_L3_MEM_TAG_WIDTH : (`L2_MEM_TAG_WIDTH + `CLOG2(`L3_NUM_REQS)))

////////////////////////// Tcache Configurable Knobs //////////////////////////

`define TCACHE_ID               $sformatf("core%0d-tcache", CORE_ID)

// Word size in bytes
`define TCACHE_WORD_SIZE        4
`define TCACHE_ADDR_WIDTH       (32-`CLOG2(`TCACHE_WORD_SIZE))

// Block size in bytes
`define TCACHE_LINE_SIZE        `L1_CACHE_LINE_SIZE

// Response tag select bits       
`define TCACHE_TAG_SEL_BITS      2

// Core request tag bits
`define TCACHE_TAG_WIDTH        (`UUID_BITS + `TCACHE_TAG_SEL_BITS)

// Input request size
`define TCACHE_NUM_REQS         `NUM_THREADS

// Memory request tag bits
`define TCACHE_MEM_TAG_WIDTH    (`CLOG2(`TCACHE_NUM_BANKS) + `CLOG2(`TCACHE_MSHR_SIZE))

////////////////////////// Rcache Configurable Knobs //////////////////////////

`define RCACHE_ID               $sformatf("cluster%0d-rcache", CLUSTER_ID)

// Word size in bytes
`define RCACHE_WORD_SIZE        4
`define RCACHE_ADDR_WIDTH       (32-`CLOG2(`RCACHE_WORD_SIZE))

// Block size in bytes
`define RCACHE_LINE_SIZE        `L2_CACHE_LINE_SIZE

// Response tag select bits       
`define RCACHE_TAG_SEL_BITS     2

// Core request tag bits
`define RCACHE_TAG_WIDTH        `RCACHE_TAG_SEL_BITS

// Input request size
`define RCACHE_NUM_REQS         1

// Memory request data bits
`define RCACHE_MEM_DATA_WIDTH   (`RCACHE_LINE_SIZE * 8)

// Memory request address bits
`define RCACHE_MEM_ADDR_WIDTH   (32 - `CLOG2(`RCACHE_LINE_SIZE))

// Memory request tag bits
`define RCACHE_MEM_TAG_WIDTH    (`CLOG2(`RCACHE_NUM_BANKS) + `CLOG2(`RCACHE_MSHR_SIZE))

////////////////////////// Ocache Configurable Knobs //////////////////////////

`define OCACHE_ID               $sformatf("cluster%0d-ocache", CLUSTER_ID)

// Word size in bytes
`define OCACHE_WORD_SIZE        4
`define OCACHE_ADDR_WIDTH       (32-`CLOG2(`OCACHE_WORD_SIZE))

// Block size in bytes
`define OCACHE_LINE_SIZE        `L2_CACHE_LINE_SIZE

// Input request size
`define OCACHE_NUM_REQS         (2*`NUM_THREADS)

// Response tag select bits       
`define OCACHE_TAG_SEL_BITS     `CLOG2(`ROP_MEM_QUEUE_SIZE)

// Core request tag bits
`define OCACHE_TAG_WIDTH        `OCACHE_TAG_SEL_BITS

// Memory request data bits
`define OCACHE_MEM_DATA_WIDTH   (`OCACHE_LINE_SIZE * 8)

// Memory request address bits
`define OCACHE_MEM_ADDR_WIDTH   (32 - `CLOG2(`OCACHE_LINE_SIZE))

// Memory request tag bits
`define OCACHE_MEM_TAG_WIDTH    (`CLOG2(`OCACHE_NUM_BANKS) + `CLOG2(`OCACHE_MSHR_SIZE))

///////////////////////////////////////////////////////////////////////////////

`define VX_MEM_BYTEEN_WIDTH     `L3_MEM_BYTEEN_WIDTH   
`define VX_MEM_ADDR_WIDTH       `L3_MEM_ADDR_WIDTH
`define VX_MEM_DATA_WIDTH       `L3_MEM_DATA_WIDTH
`define VX_MEM_TAG_WIDTH        `L3_MEM_TAG_WIDTH
`define VX_DCR_ADDR_WIDTH       `DCR_ADDR_BITS
`define VX_DCR_DATA_WIDTH       32

`define TO_FULL_ADDR(x)         {x, (32-$bits(x))'(0)}

///////////////////////////////////////////////////////////////////////////////

`define ASSIGN_VX_MEM_REQ_IF(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.rw     = src.rw;     \
    assign dst.byteen = src.byteen; \
    assign dst.addr   = src.addr;   \
    assign dst.data   = src.data;   \
    assign dst.tag    = src.tag;    \
    assign src.ready  = dst.ready

`define ASSIGN_VX_MEM_RSP_IF(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.data   = src.data;   \
    assign dst.tag    = src.tag;    \
    assign src.ready  = dst.ready

`define ASSIGN_VX_MEM_REQ_IF_XTAG(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.rw     = src.rw;     \
    assign dst.byteen = src.byteen; \
    assign dst.addr   = src.addr;   \
    assign dst.data   = src.data;   \
    assign src.ready  = dst.ready

`define ASSIGN_VX_MEM_RSP_IF_XTAG(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.data   = src.data;   \
    assign src.ready  = dst.ready

`define ASSIGN_VX_CACHE_REQ_IF(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.rw     = src.rw;     \
    assign dst.byteen = src.byteen; \
    assign dst.addr   = src.addr;   \
    assign dst.data   = src.data;   \
    assign dst.tag    = src.tag;    \
    assign src.ready  = dst.ready

`define ASSIGN_VX_CACHE_RSP_IF(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.data   = src.data;   \
    assign dst.tag    = src.tag;    \
    assign src.ready  = dst.ready

`define ASSIGN_VX_CACHE_REQ_IF_XTAG(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.rw     = src.rw;     \
    assign dst.byteen = src.byteen; \
    assign dst.addr   = src.addr;   \
    assign dst.data   = src.data;   \
    assign src.ready  = dst.ready

`define ASSIGN_VX_CACHE_RSP_IF_XTAG(dst, src) \
    assign dst.valid  = src.valid;  \
    assign dst.data   = src.data;   \
    assign src.ready  = dst.ready

`define CACHE_REQ_TO_MEM(dst, src, i) \
    assign dst[i].valid = src.valid[i]; \
    assign dst[i].rw = src.rw[i]; \
    assign dst[i].byteen = src.byteen[i]; \
    assign dst[i].addr = src.addr[i]; \
    assign dst[i].data = src.data[i]; \
    assign dst[i].tag = src.tag[i]; \
    assign src.ready[i] = dst[i].ready

`define CACHE_RSP_FROM_MEM(dst, src, i) \
    assign dst.valid[i] = src[i].valid; \
    assign dst.data[i] = src[i].data; \
    assign dst.tag[i] = src[i].tag; \
    assign src[i].ready = dst.ready[i]

///////////////////////////////////////////////////////////////////////////////

`include "VX_fpu_types.vh"
`include "VX_gpu_types.vh"

`endif
