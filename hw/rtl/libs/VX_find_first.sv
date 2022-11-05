`include "VX_platform.vh"

`TRACING_OFF
module VX_find_first #(
    parameter N       = 1,
    parameter DATAW   = 1,
    parameter REVERSE = 0
    
) (
    input  wire [N-1:0][DATAW-1:0] data_in,
    input  wire [N-1:0]            valid_in,    
    output wire [DATAW-1:0]        data_out,
    output wire                    valid_out
);
    localparam LOGN = $clog2(N);
    localparam TL   = (1 << LOGN) - 1;
    localparam TN   = (1 << (LOGN+1)) - 1;

`IGNORE_WARNINGS_BEGIN
    wire [TN-1:0] s_n;
    wire [TN-1:0][DATAW-1:0] d_n;
`IGNORE_WARNINGS_END

    for (genvar i = 0; i < N; ++i) begin
        assign s_n[TL+i] = REVERSE ? valid_in[N-1-i] : valid_in[i];
        assign d_n[TL+i] = REVERSE ? data_in[N-1-i] : data_in[i];
    end
    
    for (genvar i = TL+N; i < TN; ++i) begin
        assign s_n[i] = 0;
        assign d_n[i] = '0;
    end

    for (genvar j = 0; j < LOGN; ++j) begin
        for (genvar i = 0; i < (2**j); ++i) begin
            assign s_n[2**j-1+i] = s_n[2**(j+1)-1+i*2] | s_n[2**(j+1)-1+i*2+1];
            assign d_n[2**j-1+i] = s_n[2**(j+1)-1+i*2] ? d_n[2**(j+1)-1+i*2] : d_n[2**(j+1)-1+i*2+1];
        end
    end     
        
    assign valid_out = s_n[0];
    assign data_out  = d_n[0];
  
endmodule
`TRACING_ON
