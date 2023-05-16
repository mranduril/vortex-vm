`include "VX_define.vh"

interface VX_commit_if();

    wire                        valid;
    wire [`UP(`UUID_BITS)-1:0]  uuid;
    wire [`UP(`NW_BITS)-1:0]    wid;
    wire [`NUM_THREADS-1:0]     tmask;    
    wire [`XLEN-1:0]            PC;
    wire [`NUM_THREADS-1:0][`XLEN-1:0] data;
    wire [`NR_BITS-1:0]         rd;
    wire                        wb;
    wire                        eop;
    wire                        ready;

    modport master (
        output valid,
        output uuid,
        output wid,
        output tmask,
        output PC,
        output data,
        output rd,
        output wb,
        output eop,
        input  ready
    );

    modport slave (
        input  valid,
        input  uuid,
        input  wid,
        input  tmask,
        input  PC,
        input  data,
        input  rd,
        input  wb,
        input  eop,
        output ready
    );

endinterface
