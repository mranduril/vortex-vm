`include "VX_rop_define.vh"

module VX_rop_ds #(
    parameter DEPTH_TEST = 1,
    parameter STENCIL_TEST = 1,
    parameter CLUSTER_ID = 0,
    parameter NUM_LANES  = 4
) (
    input wire clk,
    input wire reset,

    // Depth Test
    input wire [`ROP_DEPTH_FUNC_BITS] depth_func,
    input wire [23:0] depth_ref,
    input wire [23:0] depth_val,
    input wire [3:0]  depth_mask,

    output wire [23:0] depth_out,

    // Stencil Test
    input wire [7:0]                        stencil_val;
    input wire [`ROP_STENCIL_FACE_BITS-1:0] stencil_face;

    input wire [`ROP_DEPTH_FUNC_BITS-1:0] stencil_front_func;    
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_front_zpass;
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_front_zfail;
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_front_fail;
    input wire [7:0]                      stencil_front_mask;
    input wire [7:0]                      stencil_front_ref;
    input wire [`ROP_DEPTH_FUNC_BITS-1:0] stencil_back_func;    
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_back_zpass;
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_back_zfail;
    input wire [`ROP_STENCIL_OP_BITS-1:0] stencil_back_fail;
    input wire [7:0]                      stencil_back_mask;
    input wire [7:0]                      stencil_back_ref;

    output wire [7:0]  stencil_out
    output wire [31:0] mask_out

);

    wire [`ROP_DEPTH_FUNC_BITS-1:0] stencil_func;
    reg  [`ROP_STENCIL_OP_BITS-1:0] stencil_op;
    wire [7:0] stencil_mask;
    wire [7:0] stencil_ref;

    wire [7:0] stencil_ref_m;
    wire [7:0] stencil_val_m;

    wire [23:0] depth_result;
    wire [7:0]  stencil_result;
    reg  [31:0] write_mask = stencil_mask << 24;

    wire dpass;
    wire spass;

    ///////////////////////////////////////////////////////////////

    // Depth Test

    VX_rop_compare #(
        .DATAW (24)
    ) depth_compare (
        .func   (depth_func),
        .a      (depth_ref),
        .b      (depth_val),
        .result (dpass)
    );

    always @(*) begin
        if (!DEPTH_TEST) 
            depth_result = depth_val;
        else begin
            if (dpass & (| depth_mask))
                depth_result = depth_ref;
            else    
                depth_result = depth_val;
        end
    end

    ///////////////////////////////////////////////////////////////

    // Stencil Test

    assign stencil_func = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_func : stencil_back_func;
    assign stencil_mask = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_mask : stencil_back_mask;
    assign stencil_ref  = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_ref : stencil_back_ref;

    assign stencil_ref_m = stencil_ref & stencil_mask;
    assign stencil_val_m = stencil_val & stencil_mask;

    VX_rop_compare #(
        .DATAW (8)
    ) stencil_compare (
        .func   (stencil_func),
        .a      (stencil_ref_m),
        .b      (stencil_val_m),
        .result (spass)
    );

    always @(*) begin
        if (!STENCIL_TEST)
            stencil_out = stencil_val;
        else begin 
            if (spass) begin
                if (dpass) begin
                    if (| depth_mask)
                        write_mask = write_mask | 24'hFFFFFF;
                    stencil_op = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_zpass : stencil_back_zpass;
                end else
                    stencil_op = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_zfail : stencil_back_zfail;
            end else
                stencil_op = (stencil_face == `ROP_STENCIL_FACE_FRONT)? stencil_front_fail : stencil_back_fail;
        end
    end

    VX_rop_stencil_op #(
        .DATAW (8)
    ) stencil_op_ (
        .stencil_op     (stencil_op),
        .stencil_ref    (stencil_ref),
        .stencil_val    (stencil_val),
        .stencil_result (stencil_result)
    );

    ///////////////////////////////////////////////////////////////

    VX_pipe_register #(
        .DATAW	(8 + 24 + 32),
        .RESETW (1)
    ) pipe_reg (
        .clk      (clk),
        .reset    (reset),
        .enable   (1'b1),
        .data_in  ({depth_result, stencil_result, write_mask}),
        .data_out ({depth_out,    stencil_out,    mask_out})
    );

    ///////////////////////////////////////////////////////////////

endmodule