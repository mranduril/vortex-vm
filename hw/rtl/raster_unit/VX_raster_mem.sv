`include "VX_raster_define.vh"

// Memory interface for the rasterization unit.
// Performs the following:
//  1. Break the request in tile and primitive fetch requests
//  2. Forms an FSM to keep a track of the return value types
//  3. Sotres the data in a reservation station
//  4. Dispatches the packets from the reservation station to the rasterslice
//      depending on the ready signal from the slices

module VX_raster_mem #(  
    parameter RASTER_SLICE_NUM  = 4,
    parameter RASTER_TILE_SIZE  = 16,
    parameter RASTER_TILE_BITS  = $clog2(RASTER_TILE_SIZE),
    parameter RASTER_RS_SIZE    = 8, // Reservation station size
    parameter RASTER_SLICE_BITS = `LOG2UP(RASTER_SLICE_NUM)
) (
    // Standard inputs
    input logic         clk,
    input logic         reset,

    // To indicate valid inputs provided
    input logic         input_valid,
    // Memory information
    input logic [`RASTER_DCR_DATA_BITS-1:0]         num_tiles,
    input logic [`RASTER_DCR_DATA_BITS-1:0]         tbuf_baseaddr,
    input logic [`RASTER_DCR_DATA_BITS-1:0]         pbuf_baseaddr,
    input logic [`RASTER_DCR_DATA_BITS-1:0]         pbuf_stride,

    // Raster slice interactions
    input logic [RASTER_SLICE_NUM-1:0]                          raster_slice_ready,
    output logic [`RASTER_TILE_DATA_BITS-1:0]                   out_x_loc, out_y_loc,
    output logic [`RASTER_PRIMITIVE_DATA_BITS-1:0]              out_edges[2:0][2:0],
    output logic [RASTER_SLICE_BITS-1:0]                        out_slice_index,

    // Status signals
    output logic                                    ready, out_valid
);

    logic [2:0] location_identifier;

    // Bit tag identifier for type of memory request
    // Request #1: Tile and primitive id fetch => bit = 01
    // Request #2: Primitive id fetch => bit = 10
    // Request #3: Primitive data fetch => bit = 11
    localparam TILE_FETCH      = 2'b01;
    localparam PRIM_ID_FETCH   = 2'b10;
    localparam PRIM_DATA_FETCH = 2'b11;
    localparam TILE_FETCH_MASK = 9'(2'b11);
    localparam PRIM_ID_FETCH_MASK = 9'(1'b1);
    localparam PRIM_DATA_FETCH_MASK = {9{1'b1}};

    // Temp storage to cycle through all primitives and tiles
    logic [`RASTER_DCR_DATA_BITS-1:0] temp_tbuf_addr, temp_pbuf_addr, temp_num_prims, temp_num_tiles, temp_pbuf_stride;
    logic [`RASTER_TILE_DATA_BITS-1:0] temp_x_loc, temp_y_loc;
    logic [`RASTER_DCR_DATA_BITS-1:0] temp_prim_count, temp_tile_count;

    // Holds x_loc, y_loc, edge_func_val, edges -> extents are calculated on the fly
    localparam RASTER_RS_DATA_WIDTH = 2*`RASTER_TILE_DATA_BITS + 3*3*`RASTER_PRIMITIVE_DATA_BITS;
    localparam RASTER_RS_INDEX_BITS = `LOG2UP(RASTER_RS_SIZE);

    // Reservation station
    logic [RASTER_RS_DATA_WIDTH-1:0]    raster_rs[RASTER_RS_SIZE-1:0];
    logic [RASTER_RS_SIZE-1:0]          raster_rs_valid;
    logic [RASTER_RS_SIZE-1:0]          raster_rs_empty;

    logic [RASTER_RS_INDEX_BITS-1:0]    raster_rs_empty_index, raster_rs_index;

    // Status signals
    logic valid_raster_index, valid_rs_index, valid_rs_empty_index;
    logic fetch_fsm_complete;

    // Memory interactions
    logic mem_req_valid, mem_req_ready, mem_rsp_valid;
    logic [8:0][`RASTER_DCR_DATA_BITS-1:0]       mem_req_addr;
    logic [8:0][`RASTER_PRIMITIVE_DATA_BITS-1:0] mem_rsp_data;
    localparam TAG_MAX_BIT_INDEX = RASTER_RS_INDEX_BITS-1 + 2;
    logic [TAG_MAX_BIT_INDEX:0]                  mem_rsp_tag; // size increased by 1 bit to account for the mem_tag_type
    logic [1:0]                                  mem_tag_type;
    logic [8:0]                                  mem_req_mask;

    // Stall signal
    //  -> assert when any entry in the RS is empty
    assign ready = |raster_rs_empty & mem_req_ready & fetch_fsm_complete;

    always @(posedge clk) begin
        // Reset condition
        if (reset) begin
            mem_req_valid <= 0;
            out_valid <= 0;
            temp_tbuf_addr <= 0;
            temp_pbuf_addr <= 0;
            temp_num_prims <= 0;
            temp_num_tiles <= 0;
            temp_prim_count <= 0;
            fetch_fsm_complete <= 1;
            for (int i = 0; i < RASTER_RS_SIZE; ++i) begin
                raster_rs_valid[i] <= 0;
                raster_rs_empty[i] <= 1;
            end
        end
        else if (mem_req_ready) begin
            // On new input -> set the temp state values
            if (ready && input_valid && valid_rs_empty_index) begin
                temp_pbuf_addr <= pbuf_baseaddr;
                temp_num_tiles <= num_tiles;
                temp_pbuf_stride <= pbuf_stride;
                // Launch memory request to get the tbuf and pbuf
                mem_tag_type <= TILE_FETCH;
                mem_req_mask <= TILE_FETCH_MASK;
                mem_req_valid <= 1;
                // Fetch the first primitive as well
                mem_req_addr[0] <= tbuf_baseaddr;
                mem_req_addr[1] <= tbuf_baseaddr + 4;
                // To indicate the address for next fetch
                temp_tbuf_addr <= tbuf_baseaddr + 4 + 4;

                // Reset the counter
                temp_prim_count <= 0;
                temp_tile_count <= 0;
                fetch_fsm_complete <= 0;
            end
            // If it gets a valid memory response
            if (mem_rsp_valid && mem_req_ready) begin
                // Check if all tiles have been generated
                if (temp_tile_count == temp_num_tiles) begin
                    fetch_fsm_complete <= 1;
                end
                // If not generate the tiles and primitives
                else begin
                    // Check the reponse tag type
                    if (mem_rsp_tag[1:0] == TILE_FETCH) begin
                        location_identifier <= 1;
                        // returned value is tile data
                        temp_x_loc <= `RASTER_TILE_DATA_BITS'((mem_rsp_data[0] & {{`RASTER_TILE_DATA_BITS{1'b0}}, {`RASTER_TILE_DATA_BITS{1'b1}}}) << RASTER_TILE_BITS);
                        temp_y_loc <= `RASTER_TILE_DATA_BITS'((mem_rsp_data[0] >> (`RASTER_TILE_DATA_BITS)) << RASTER_TILE_BITS);
                        temp_num_prims <= mem_rsp_data[1];
                        // Launch the pid fetch
                        mem_tag_type <= PRIM_ID_FETCH;
                        mem_req_valid <= 1;
                        mem_req_mask <= PRIM_ID_FETCH_MASK;
                        // Fetch the primitive index
                        mem_req_addr[0] <= temp_tbuf_addr;
                        temp_tbuf_addr <= temp_tbuf_addr + 4;
                    end
                    else if (mem_rsp_tag[1:0] == PRIM_ID_FETCH) begin
                        location_identifier <= 2;
                        //integer pbuf_addr = temp_pbuf_addr + mem_rsp_data[0]*temp_pbuf_stride;
                        // Launch next request based on pid
                        mem_tag_type <= PRIM_DATA_FETCH;
                        mem_req_valid <= 1;
                        mem_req_mask <= PRIM_DATA_FETCH_MASK;
                        for (int i = 0; i < 9; ++i) begin
                            mem_req_addr[i] <= temp_pbuf_addr + mem_rsp_data[0] * temp_pbuf_stride + 4*i;
                        end
                    end
                    else if (mem_rsp_tag[1:0] == PRIM_DATA_FETCH) begin
                        location_identifier <= 3;
                        // Insert data into RS
                        raster_rs[raster_rs_empty_index] <= {temp_x_loc, temp_y_loc,
                            mem_rsp_data[0], mem_rsp_data[1], mem_rsp_data[2],
                            mem_rsp_data[3], mem_rsp_data[4], mem_rsp_data[5],
                            mem_rsp_data[6], mem_rsp_data[7], mem_rsp_data[8]
                        };
                        raster_rs_empty[raster_rs_empty_index] <= 0;
                        raster_rs_valid[raster_rs_empty_index] <= 1;

                        // Incrememnt the prim and tile count register
                        if (temp_prim_count + 1 == temp_num_prims) begin
                            location_identifier <= 4;
                            // => Last primitive fetched
                            temp_prim_count <= 0;
                            temp_tile_count <= temp_tile_count + 1;
                            // Check if this was last tile
                            if (temp_tile_count + 1 == temp_num_tiles) begin
                                // => this was last tile
                                temp_tile_count <= temp_tile_count + 1;
                            end
                            else begin
                                location_identifier <= 5;
                                // Fetch the next tile
                                // Launch memory request to get the tbuf and pbuf
                                mem_tag_type <= TILE_FETCH;
                                mem_req_valid <= 1;
                                mem_req_mask <= TILE_FETCH_MASK;
                                // Fetch the first primitive as well
                                mem_req_addr[0] <= temp_tbuf_addr;
                                mem_req_addr[1] <= temp_tbuf_addr + 4;
                                // To indicate the address for next fetch
                                temp_tbuf_addr <= temp_tbuf_addr + 4 + 4;
                            end
                        end
                        else begin
                            location_identifier <= 6;
                            temp_prim_count <= temp_prim_count + 1;
                            // Launch the request to get next primitive id // Launch the pid fetch
                            mem_tag_type <= PRIM_ID_FETCH;
                            mem_req_valid <= 1;
                            mem_req_mask <= PRIM_ID_FETCH_MASK;
                            // Fetch the primitive index
                            mem_req_addr[0] <= temp_tbuf_addr;
                            temp_tbuf_addr <= temp_tbuf_addr + 4;
                        end
                    end
                    else begin
                        `ASSERT(0, ("Incorrect tag returned"));
                        mem_req_valid <= 0;
                    end
                end
            end

            // Launch any valid packet
            // When any raster slice is ready
            if (valid_raster_index && valid_rs_index) begin
                {out_x_loc, out_y_loc,
                out_edges[0][0], out_edges[0][1], out_edges[0][2],
                out_edges[1][0], out_edges[1][1], out_edges[1][2],
                out_edges[2][0], out_edges[2][1], out_edges[2][2]} <= raster_rs[raster_rs_index];
                raster_rs_valid[raster_rs_index] <= 0;
                raster_rs_empty[raster_rs_index] <= 1;
                out_valid <= 1;
            end
            else begin
                out_valid <= 0;
            end
        end
    end

    // Priority encoder to select a free index in the RS
    VX_priority_encoder #( 
        .N(RASTER_RS_SIZE)
    ) raster_empty_rs (
        .data_in(raster_rs_empty),  
        `UNUSED_PIN (onehot),
        .index(raster_rs_empty_index),
        .valid_out(valid_rs_empty_index)
    );

    // Priority encoder to select the valid entry from RS to dispatch
    VX_priority_encoder #( 
        .N(RASTER_RS_SIZE)
    ) raster_request_rs (
        .data_in(raster_rs_valid),  
        `UNUSED_PIN (onehot),
        .index(raster_rs_index),
        .valid_out(valid_rs_index)
    );

    VX_priority_encoder #(
        .N(RASTER_SLICE_NUM)
    ) raster_ready_select (
        .data_in(raster_slice_ready),
        `UNUSED_PIN (onehot),
        .index(out_slice_index),
        .valid_out(valid_raster_index)
    );

    // Memory streamer
    VX_mem_streamer #(
        .NUM_REQS(9), // 3 edges and 3 coeffs in each edge
        .ADDRW(`RASTER_DCR_DATA_BITS),
        .DATAW(`RASTER_PRIMITIVE_DATA_BITS),
        .TAGW(RASTER_RS_INDEX_BITS + 2) // the top bit will denote type of request
    ) raster_mem_streamer (
        .clk(clk),
        .reset(reset),

        .req_valid(mem_req_valid && mem_req_ready && |raster_rs_empty), // NOTE: This should ensure stalls
        .req_rw(0),
        .req_mask(mem_req_mask),
        `UNUSED_PIN (req_byteen),   /// TODO: USE THIS PIN
        .req_addr(mem_req_addr),
        `UNUSED_PIN (req_data),
        .req_tag({raster_rs_empty_index, mem_tag_type}), //tag type appended to tag
        .req_ready(mem_req_ready),
        
        // Output response
        .rsp_valid(mem_rsp_valid),
        `UNUSED_PIN (rsp_mask),
        .rsp_data(mem_rsp_data),
        .rsp_tag(mem_rsp_tag),
        .rsp_ready(1),

        `UNUSED_PIN (mem_req_valid),
        `UNUSED_PIN (mem_req_rw),
        `UNUSED_PIN (mem_req_byteen),
        `UNUSED_PIN (mem_req_addr),
        `UNUSED_PIN (mem_req_data),
        `UNUSED_PIN (mem_req_tag),
        `UNUSED_PIN (mem_req_ready),
        `UNUSED_PIN (mem_rsp_valid),
        `UNUSED_PIN (mem_rsp_mask),
        `UNUSED_PIN (mem_rsp_data),
        `UNUSED_PIN (mem_rsp_tag),
        `UNUSED_PIN (mem_rsp_ready)
    );

endmodule