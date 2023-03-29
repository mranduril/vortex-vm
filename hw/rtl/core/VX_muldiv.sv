`include "VX_define.vh"

module VX_muldiv (
    input wire clk,
    input wire reset,
    
    // Inputs    
    input wire [`INST_MUL_BITS-1:0]     alu_op,
    input wire [`UP(`UUID_BITS)-1:0]    uuid_in,
    input wire [`UP(`NW_BITS)-1:0]      wid_in,
    input wire [`NUM_THREADS-1:0]       tmask_in,
    input wire [31:0]                   PC_in,
    input wire [`NR_BITS-1:0]           rd_in,
    input wire                          wb_in,
    input wire [`NUM_THREADS-1:0][`XLEN-1:0] alu_in1, 
    input wire [`NUM_THREADS-1:0][`XLEN-1:0] alu_in2,

    // Outputs
    output wire [`UP(`UUID_BITS)-1:0]    uuid_out,
    output wire [`UP(`NW_BITS)-1:0]      wid_out,
    output wire [`NUM_THREADS-1:0]       tmask_out,
    output wire [31:0]                   PC_out,
    output wire [`NR_BITS-1:0]           rd_out,
    output wire                          wb_out,
    output wire [`NUM_THREADS-1:0][`XLEN-1:0] data_out,

    // handshake
    input wire  valid_in,
    output wire ready_in,
    output wire valid_out,
    input wire  ready_out
);
    localparam UUID_WIDTH = `UP(`UUID_BITS);
    localparam NW_WIDTH   = `UP(`NW_BITS); 

    wire is_div_op = `INST_MUL_IS_DIV(alu_op);

    wire [`NUM_THREADS-1:0][`XLEN-1:0] mul_result;
    wire [UUID_WIDTH-1:0] mul_uuid_out;
    wire [NW_WIDTH-1:0] mul_wid_out;
    wire [`NUM_THREADS-1:0] mul_tmask_out;
    wire [31:0] mul_PC_out;
    wire [`NR_BITS-1:0] mul_rd_out;
    wire mul_wb_out;

    wire stall_out;

    wire mul_valid_out;
    wire mul_valid_in = valid_in && !is_div_op;    
    wire mul_ready_in = ~stall_out || ~mul_valid_out;

    wire is_mulh_in      = (alu_op != `INST_MUL_MUL) && (alu_op != `INST_MUL_MULW);
    wire is_signed_mul_a = (alu_op != `INST_MUL_MULHU);
    wire is_signed_mul_b = (alu_op != `INST_MUL_MULHU && alu_op != `INST_MUL_MULHSU);
    wire is_mulw         = (alu_op == `INST_MUL_MULW);

`ifdef IMUL_DPI

    wire [`NUM_THREADS-1:0][`XLEN-1:0] mul_result_tmp;  

    wire mul_fire_in = mul_valid_in && mul_ready_in;

    for (genvar i = 0; i < `NUM_THREADS; ++i) begin
        wire [`XLEN-1:0] mul_resultl, mul_resulth;
        wire [`XLEN-1:0] alu_in1_dpi = is_mulw ? (alu_in1[i] & `XLEN'hFFFFFFFF) : alu_in1[i]; 
        wire [`XLEN-1:0] alu_in2_dpi = is_mulw ? (alu_in2[i] & `XLEN'hFFFFFFFF) : alu_in2[i]; 
        always @(*) begin   
            dpi_imul (mul_fire_in, alu_in1_dpi, alu_in2_dpi, is_signed_mul_a, is_signed_mul_b, mul_resultl, mul_resulth);
        end
        assign mul_result_tmp[i] = is_mulh_in ? mul_resulth : (is_mulw ? `XLEN'($signed(mul_resultl[31:0])) : mul_resultl);
    end

    VX_shift_register #(
        .DATAW  (1 + UUID_WIDTH + NW_WIDTH + `NUM_THREADS + 32 + `NR_BITS + 1 + (`NUM_THREADS * `XLEN)),
        .DEPTH  (`LATENCY_IMUL),
        .RESETW (1)
    ) mul_shift_reg (
        .clk(clk),
        .reset    (reset),
        .enable   (mul_ready_in),
        .data_in  ({mul_valid_in,  uuid_in,      wid_in,      tmask_in,       PC_in,      rd_in,      wb_in,      mul_result_tmp}),
        .data_out ({mul_valid_out, mul_uuid_out, mul_wid_out, mul_tmask_out,  mul_PC_out, mul_rd_out, mul_wb_out, mul_result})
    );

`else      
    
    wire is_mulh_out;
    //TODO handle mulw when not using DPI
    for (genvar i = 0; i < `NUM_THREADS; ++i) begin
        wire [`XLEN:0] mul_in1 = {is_signed_mul_a && alu_in1[i][`XLEN-1], alu_in1[i]};
        wire [`XLEN:0] mul_in2 = {is_signed_mul_b && alu_in2[i][`XLEN-1], alu_in2[i]};
    `IGNORE_UNUSED_BEGIN
        wire [2*(`XLEN)+1:0] mul_result_tmp;
    `IGNORE_UNUSED_END

        VX_multiplier #(
            .A_WIDTH (`XLEN+1),
            .B_WIDTH (`XLEN+1),
            .R_WIDTH (2*(`XLEN)+1),
            .SIGNED  (1),
            .LATENCY (`LATENCY_IMUL)
        ) multiplier (
            .clk    (clk),
            .enable (mul_ready_in),
            .dataa  (mul_in1),
            .datab  (mul_in2),
            .result (mul_result_tmp)
        );

        assign mul_result[i] = is_mulh_out ? mul_result_tmp[2*(`XLEN)-1:`XLEN] : mul_result_tmp[`XLEN-1:0];
    end

    VX_shift_register #(
        .DATAW  (1 + UUID_WIDTH + NW_WIDTH + `NUM_THREADS + 32 + `NR_BITS + 1 + 1),
        .DEPTH  (`LATENCY_IMUL),
        .RESETW (1)
    ) mul_shift_reg (
        .clk(clk),
        .reset    (reset),
        .enable   (mul_ready_in),
        .data_in  ({mul_valid_in,  uuid_in,      wid_in,      tmask_in,       PC_in,      rd_in,      wb_in,      is_mulh_in}),
        .data_out ({mul_valid_out, mul_uuid_out, mul_wid_out, mul_tmask_out,  mul_PC_out, mul_rd_out, mul_wb_out, is_mulh_out})
    );

`endif

    ///////////////////////////////////////////////////////////////////////////

    wire [`NUM_THREADS-1:0][`XLEN-1:0] div_result;
    wire [UUID_WIDTH-1:0] div_uuid_out;
    wire [NW_WIDTH-1:0] div_wid_out;
    wire [`NUM_THREADS-1:0] div_tmask_out;
    wire [31:0] div_PC_out;
    wire [`NR_BITS-1:0] div_rd_out;
    wire div_wb_out;

    wire is_rem_op_in  = (alu_op == `INST_MUL_REM) || (alu_op == `INST_MUL_REMU);
    wire is_signed_div = (alu_op == `INST_MUL_DIV) || (alu_op == `INST_MUL_REM) || (alu_op == `INST_MUL_DIVW);     
    wire div_valid_in  = valid_in && is_div_op; 
    wire div_ready_out = ~stall_out && ~mul_valid_out; // arbitration prioritizes MUL  
    wire div_ready_in;
    wire div_valid_out;

    wire is_divw         = (alu_op == `INST_MUL_DIVW);
    wire is_divuw        = (alu_op == `INST_MUL_DIVUW);

`ifdef IDIV_DPI    

    wire [`NUM_THREADS-1:0][`XLEN-1:0] div_result_tmp;

    wire div_fire_in = div_valid_in && div_ready_in;
    for (genvar i = 0; i < `NUM_THREADS; ++i) begin
        wire [`XLEN-1:0] div_quotient, div_remainder;
        wire [`XLEN-1:0] alu_in1_dpi = is_divuw ? (alu_in1[i] & `XLEN'hFFFFFFFF) : (is_divw ? `XLEN'($signed(alu_in1[i][31:0])): alu_in1[i]);
        wire [`XLEN-1:0] alu_in2_dpi = is_divuw ? (alu_in2[i] & `XLEN'hFFFFFFFF) : (is_divw ? `XLEN'($signed(alu_in2[i][31:0])): alu_in2[i]);
        always @(*) begin  
            dpi_idiv (div_fire_in, alu_in1_dpi, alu_in2_dpi, is_signed_div, div_quotient, div_remainder);
        end
        wire [`XLEN-1:0] div_quotient_out = (is_divuw | is_divw) ? `XLEN'($signed(div_quotient[31:0])) : div_quotient;
        assign div_result_tmp[i] = is_rem_op_in ? div_remainder : div_quotient_out;
    end

    VX_shift_register #(
        .DATAW  (1 + UUID_WIDTH + NW_WIDTH + `NUM_THREADS + 32 + `NR_BITS + 1 + (`NUM_THREADS * `XLEN)),
        .DEPTH  (`LATENCY_IMUL),
        .RESETW (1)
    ) div_shift_reg (
        .clk(clk),
        .reset    (reset),
        .enable   (div_ready_in),
        .data_in  ({div_valid_in,  uuid_in,      wid_in,      tmask_in,       PC_in,      rd_in,      wb_in,      div_result_tmp}),
        .data_out ({div_valid_out, div_uuid_out, div_wid_out, div_tmask_out,  div_PC_out, div_rd_out, div_wb_out, div_result})
    );

    assign div_ready_in = div_ready_out || ~div_valid_out;

`else

    wire [`NUM_THREADS-1:0][`XLEN-1:0] div_result_tmp, rem_result_tmp;
    wire is_rem_op_out;

    VX_serial_div #(
        .WIDTHN (32),
        .WIDTHD (32),
        .WIDTHQ (32),
        .WIDTHR (32),
        .LANES  (`NUM_THREADS),
        .TAGW   (UUID_WIDTH + NW_WIDTH + `NUM_THREADS + 32 + `NR_BITS + 1 + 1)
    ) divide (
        .clk       (clk),
        .reset     (reset),
        
        .valid_in  (div_valid_in),
        .ready_in  (div_ready_in),
        .tag_in    ({uuid_in, wid_in, tmask_in, PC_in, rd_in, wb_in, is_rem_op_in}),

        .ready_out (div_ready_out),
        .valid_out (div_valid_out),
        .tag_out   ({div_uuid_out, div_wid_out, div_tmask_out, div_PC_out, div_rd_out, div_wb_out, is_rem_op_out}),

        .is_signed (is_signed_div),        
        .numer     (alu_in1),
        .denom     (alu_in2),

        .quotient  (div_result_tmp),
        .remainder (rem_result_tmp)        
    );

    assign div_result = is_rem_op_out ? rem_result_tmp : div_result_tmp; 

`endif

    ///////////////////////////////////////////////////////////////////////////

    wire                    rsp_valid = mul_valid_out || div_valid_out;  
    wire [UUID_WIDTH-1:0]   rsp_uuid = mul_valid_out ? mul_uuid_out : div_uuid_out;
    wire [NW_WIDTH-1:0]     rsp_wid  = mul_valid_out ? mul_wid_out : div_wid_out;
    wire [`NUM_THREADS-1:0] rsp_tmask = mul_valid_out ? mul_tmask_out : div_tmask_out;
    wire [31:0]             rsp_PC    = mul_valid_out ? mul_PC_out : div_PC_out;
    wire [`NR_BITS-1:0]     rsp_rd    = mul_valid_out ? mul_rd_out : div_rd_out;
    wire                    rsp_wb    = mul_valid_out ? mul_wb_out : div_wb_out;
    wire [`NUM_THREADS-1:0][`XLEN-1:0] rsp_data = mul_valid_out ? mul_result : div_result;

    assign stall_out = ~ready_out && valid_out;

    VX_pipe_register #(
        .DATAW  (1 + UUID_WIDTH + NW_WIDTH + `NUM_THREADS + 32 + `NR_BITS + 1 + (`NUM_THREADS * `XLEN)),
        .RESETW (1)
    ) pipe_reg (
        .clk      (clk),
        .reset    (reset),
        .enable   (~stall_out),
        .data_in  ({rsp_valid, rsp_uuid, rsp_wid, rsp_tmask, rsp_PC, rsp_rd, rsp_wb, rsp_data}),
        .data_out ({valid_out, uuid_out, wid_out, tmask_out, PC_out, rd_out, wb_out, data_out})
    );

    // can accept new request?
    assign ready_in = is_div_op ? div_ready_in : mul_ready_in;
    
endmodule