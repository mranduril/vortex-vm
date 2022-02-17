`ifndef VX_RASTER_TYPES
`define VX_RASTER_TYPES

`include "VX_define.vh"

package raster_types;

typedef struct packed {
    logic [31:0]    pidx_addr;
    logic [31:0]    pidx_size;
    logic [31:0]    pbuf_addr;
    logic [31:0]    pbuf_stride;
    logic [15:0]    tile_left;
    logic [15:0]    tile_top;
    logic [15:0]    tile_width;
    logic [15:0]    tile_height;
} raster_csrs_t;

endpackage

`endif