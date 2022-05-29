`include "VX_tex_define.vh"

module VX_tex_unit #(  
    parameter CORE_ID = 0
) (
    input wire  clk,
    input wire  reset,    

    // PERF
`ifdef PERF_ENABLE
    VX_tex_perf_if.master   tex_perf_if,
`endif

    // Memory interface
    VX_cache_req_if.master  cache_req_if,
    VX_cache_rsp_if.slave   cache_rsp_if,

    // Inputs
    VX_tex_dcr_if.slave     tex_dcr_if,
    VX_gpu_csr_if.slave     tex_csr_if,
    VX_tex_req_if.slave     tex_req_if,
    
    // Outputs
    VX_commit_if.master     tex_rsp_if
);

    localparam REQ_INFO_W = `NR_BITS + 1 + `NW_BITS + 32 + `UUID_BITS;
    localparam BLEND_FRAC_W = (2 * `NUM_THREADS * `TEX_BLEND_FRAC);    

    wire [`NUM_THREADS-1:0][`TEX_LOD_BITS-1:0]      mip_level;
    wire [`NUM_THREADS-1:0][`TEX_MIPOFF_BITS-1:0]   sel_mipoff;
    wire [`NUM_THREADS-1:0][1:0][`TEX_LOD_BITS-1:0] sel_logdims;

    // CSRs access

    tex_csrs_t tex_csrs;

    VX_tex_csr #(
        .CORE_ID    (CORE_ID),
        .NUM_STAGES (`TEX_STAGE_COUNT)
    ) tex_csr (
        .clk        (clk),
        .reset      (reset),

        // inputs
        .tex_csr_if (tex_csr_if),

        // outputs
        .tex_csrs   (tex_csrs)
    );

    tex_dcrs_t tex_dcrs;
    always @(posedge clk) begin
        tex_dcrs <= tex_dcr_if.data[tex_csrs.stage];
    end

    // mipmap select

    for (genvar i = 0; i < `NUM_THREADS; ++i) begin
        assign mip_level[i]      = tex_req_if.lod[i][`TEX_LOD_BITS-1:0];
        assign sel_mipoff[i]     = tex_dcrs.mipoff[mip_level[i]];
        assign sel_logdims[i][0] = tex_dcrs.logdims[0];
        assign sel_logdims[i][1] = tex_dcrs.logdims[1];
    end

    // address generation

    wire mem_req_valid;
    wire [`NUM_THREADS-1:0] mem_req_tmask;
    wire [`TEX_FILTER_BITS-1:0] mem_req_filter;
    wire [`TEX_LGSTRIDE_BITS-1:0] mem_req_lgstride;
    wire [`NUM_THREADS-1:0][1:0][`TEX_BLEND_FRAC-1:0] mem_req_blends;
    wire [`NUM_THREADS-1:0][3:0][31:0] mem_req_addr;
    wire [`NUM_THREADS-1:0][31:0] mem_req_baseaddr;
    wire [(`TEX_FORMAT_BITS + REQ_INFO_W)-1:0] mem_req_info;
    wire mem_req_ready;
                
    VX_tex_addr #(
        .CORE_ID   (CORE_ID),
        .REQ_INFOW (`TEX_FORMAT_BITS + REQ_INFO_W),
        .NUM_REQS  (`NUM_THREADS)
    ) tex_addr (
        .clk        (clk),
        .reset      (reset),

        // inputs
        .req_valid  (tex_req_if.valid),
        .req_tmask  (tex_req_if.tmask),
        .req_coords (tex_req_if.coords),
        .req_format (tex_dcrs.format),
        .req_filter (tex_dcrs.filter),
        .req_wraps  (tex_dcrs.wraps),
        .req_baseaddr(tex_dcrs.baddr),    
        .mip_level  (mip_level),
        .req_mipoff (sel_mipoff),
        .req_logdims(sel_logdims),
        .req_info   ({tex_dcrs.format, tex_req_if.rd, tex_req_if.wb, tex_req_if.wid, tex_req_if.PC, tex_req_if.uuid}),
        .req_ready  (tex_req_if.ready),

        // outputs
        .rsp_valid  (mem_req_valid), 
        .rsp_tmask  (mem_req_tmask),
        .rsp_filter (mem_req_filter), 
        .rsp_lgstride(mem_req_lgstride),
        .rsp_baseaddr(mem_req_baseaddr),
        .rsp_addr   (mem_req_addr),
        .rsp_blends (mem_req_blends),
        .rsp_info   (mem_req_info),
        .rsp_ready  (mem_req_ready)
    );

    // retrieve texel values from memory  

    wire mem_rsp_valid;
    wire [`NUM_THREADS-1:0] mem_rsp_tmask;
    wire [`NUM_THREADS-1:0][3:0][31:0] mem_rsp_data;
    wire [(BLEND_FRAC_W + `TEX_FORMAT_BITS + REQ_INFO_W)-1:0] mem_rsp_info;
    wire mem_rsp_ready;        

    VX_tex_mem #(
        .CORE_ID   (CORE_ID),
        .REQ_INFOW (BLEND_FRAC_W + `TEX_FORMAT_BITS + REQ_INFO_W),
        .NUM_REQS  (`NUM_THREADS)
    ) tex_mem (
        .clk       (clk),
        .reset     (reset),

        // memory interface
        .cache_req_if (cache_req_if),
        .cache_rsp_if (cache_rsp_if),

        // inputs
        .req_valid (mem_req_valid),
        .req_tmask (mem_req_tmask),
        .req_filter(mem_req_filter), 
        .req_lgstride(mem_req_lgstride),
        .req_baseaddr(mem_req_baseaddr),    
        .req_addr  (mem_req_addr),
        .req_info  ({mem_req_blends, mem_req_info}),
        .req_ready (mem_req_ready),

        // outputs
        .rsp_valid (mem_rsp_valid),
        .rsp_tmask (mem_rsp_tmask),
        .rsp_data  (mem_rsp_data),
        .rsp_info  (mem_rsp_info),
        .rsp_ready (mem_rsp_ready)
    );

    // apply sampler

    wire sampler_rsp_valid;
    wire [`NUM_THREADS-1:0] sampler_rsp_tmask;
    wire [`NUM_THREADS-1:0][31:0] sampler_rsp_data;
    wire [REQ_INFO_W-1:0] sampler_rsp_info;
    wire sampler_rsp_ready;

    VX_tex_sampler #(
        .CORE_ID   (CORE_ID),
        .REQ_INFOW (REQ_INFO_W),
        .NUM_REQS  (`NUM_THREADS)
    ) tex_sampler (
        .clk        (clk),
        .reset      (reset),

        // inputs
        .req_valid  (mem_rsp_valid),  
        .req_tmask  (mem_rsp_tmask),
        .req_data   (mem_rsp_data), 
        .req_blends (mem_rsp_info[(REQ_INFO_W + `TEX_FORMAT_BITS) +: BLEND_FRAC_W]),
        .req_format (mem_rsp_info[REQ_INFO_W +: `TEX_FORMAT_BITS]),
        .req_info   (mem_rsp_info[0 +: REQ_INFO_W]),
        .req_ready  (mem_rsp_ready),

        // outputs
        .rsp_valid  (sampler_rsp_valid),
        .rsp_tmask  (sampler_rsp_tmask),
        .rsp_data   (sampler_rsp_data),
        .rsp_info   (sampler_rsp_info),
        .rsp_ready  (sampler_rsp_ready)
    );

    VX_skid_buffer #(
        .DATAW (REQ_INFO_W + `NUM_THREADS * (1 + 32))
    ) rsp_sbuf (
        .clk       (clk),
        .reset     (reset),
        .valid_in  (sampler_rsp_valid),
        .ready_in  (sampler_rsp_ready),
        .data_in   ({sampler_rsp_info, sampler_rsp_tmask, sampler_rsp_data}),
        .data_out  ({tex_rsp_if.rd, tex_rsp_if.wb, tex_rsp_if.wid, tex_rsp_if.PC, tex_rsp_if.uuid, tex_rsp_if.tmask, tex_rsp_if.data}),
        .valid_out (tex_rsp_if.valid),
        .ready_out (tex_rsp_if.ready)
    );

    assign tex_rsp_if.eop = 1;

`ifdef PERF_ENABLE
    wire [$clog2(`NUM_THREADS+1)-1:0] perf_mem_req_per_cycle;
    wire [$clog2(`NUM_THREADS+1)-1:0] perf_mem_rsp_per_cycle;
    wire [$clog2(`NUM_THREADS+1)+1-1:0] perf_pending_reads_cycle;

    wire [`NUM_THREADS-1:0] perf_mem_req_per_req = cache_req_if.valid & cache_req_if.ready;
    wire [`NUM_THREADS-1:0] perf_mem_rsp_per_req = cache_rsp_if.valid & cache_rsp_if.ready;

    `POP_COUNT(perf_mem_req_per_cycle, perf_mem_req_per_req);
    `POP_COUNT(perf_mem_rsp_per_cycle, perf_mem_rsp_per_req);

    reg [`PERF_CTR_BITS-1:0] perf_pending_reads;   
    assign perf_pending_reads_cycle = perf_mem_req_per_cycle - perf_mem_rsp_per_cycle;

    always @(posedge clk) begin
        if (reset) begin
            perf_pending_reads <= 0;
        end else begin
            perf_pending_reads <= perf_pending_reads + `PERF_CTR_BITS'($signed(perf_pending_reads_cycle));
        end
    end

    reg [`PERF_CTR_BITS-1:0] perf_mem_reads;
    reg [`PERF_CTR_BITS-1:0] perf_mem_latency;

    always @(posedge clk) begin
        if (reset) begin
            perf_mem_reads   <= 0;
            perf_mem_latency <= 0;
        end else begin
            perf_mem_reads   <= perf_mem_reads + `PERF_CTR_BITS'(perf_mem_req_per_cycle);
            perf_mem_latency <= perf_mem_latency + `PERF_CTR_BITS'(perf_pending_reads);
        end
    end

    assign tex_perf_if.mem_reads   = perf_mem_reads;
    assign tex_perf_if.mem_latency = perf_pending_reads;
`endif  

`ifdef DBG_TRACE_TEX
    always @(posedge clk) begin
        if (tex_req_if.valid && tex_req_if.ready) begin
            `TRACE(1, ("%d: core%0d-tex-req: wid=%0d, PC=0x%0h, tmask=%b, stage=%0d, lod=0x%0h, u=", 
                $time, CORE_ID, tex_req_if.wid, tex_req_if.PC, tex_req_if.tmask, tex_csrs.stage, tex_req_if.lod));
            `TRACE_ARRAY1D(1, tex_req_if.coords[0], `NUM_THREADS);
            `TRACE(1, (", v="));
            `TRACE_ARRAY1D(1, tex_req_if.coords[1], `NUM_THREADS);
            `TRACE(1, (" (#%0d)\n", tex_req_if.uuid));
        end
        if (tex_rsp_if.valid && tex_rsp_if.ready) begin
             `TRACE(1, ("%d: core%0d-tex-rsp: wid=%0d, PC=0x%0h, tmask=%b, data=",
                    $time, CORE_ID, tex_rsp_if.wid, tex_rsp_if.PC, tex_rsp_if.tmask));
            `TRACE_ARRAY1D(1, tex_rsp_if.data, `NUM_THREADS);
            `TRACE(1, (" (#%0d)\n", tex_rsp_if.uuid));
        end
    end
`endif

endmodule
