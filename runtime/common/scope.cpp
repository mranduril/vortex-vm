#include "scope.h"
#include <VX_config.h>
#include <nlohmann_json.hpp>
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <vector>
#include <list>
#include <assert.h>
#include <chrono>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <unordered_set>
#include <sstream>

#define FRAME_FLUSH_SIZE 100

#define MMIO_SCOPE_READ  (AFU_IMAGE_MMIO_SCOPE_READ * 4)
#define MMIO_SCOPE_WRITE (AFU_IMAGE_MMIO_SCOPE_WRITE * 4)

#define CMD_GET_WIDTH   0
#define CMD_GET_COUNT   1
#define CMD_GET_START   2
#define CMD_GET_DATA    3
#define CMD_SET_START   4
#define CMD_SET_STOP    5

#define CHECK_ERR(_expr)    \
    do {                    \
        int err = _expr;    \
        if (err == 0)       \
            break;          \
        printf("[SCOPE] error: '%s' returned %d!\n", #_expr, err); \
        return err;         \
    } while (false)

struct tap_signal_t {
    uint32_t id;  
    std::string name;    
    uint32_t width;    
};

struct tap_t {
    uint32_t id;    
    uint32_t width;    
    uint32_t frames;    
    uint32_t cur_frame;
    uint64_t ticks;
    std::string path;
    std::vector<tap_signal_t> signals;
};

static scope_callback_t g_callback;

using json = nlohmann::json;

static uint64_t dump_clock(std::ofstream& ofs, uint64_t delta, uint64_t timestamp) {
    while (delta != 0) {
        ofs << '#' << timestamp++ << std::endl;
        ofs << "b0 0" << std::endl;
        ofs << '#' << timestamp++ << std::endl;
        ofs << "b1 0" << std::endl;
        --delta;
    }
    return timestamp;
}

static std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

static void dump_module(std::ofstream& ofs, 
                        const std::string& name,
                        std::unordered_map<std::string, std::unordered_set<std::string>>& hierarchy,
                        std::unordered_map<std::string, tap_t*>& tails,
                        int indentation) {
    std::string indent(indentation, ' ');
    ofs << indent << "$scope module " << name << " $end" << std::endl;

    auto itt = tails.find(name);
    if (itt != tails.end()) {
        for (auto& signal : itt->second->signals) {
            ofs << indent << " $var reg " << signal.width << " " << signal.id << " " << signal.name << " $end" << std::endl;                        
        }
    }

    auto ith = hierarchy.find(name);
    if (ith != hierarchy.end()) {
        for (auto& child : ith->second) {
            dump_module(ofs, child, hierarchy, tails, indentation + 1);
        }
    }

    ofs << indent << "$upscope $end" << std::endl;
}

static void dump_header(std::ofstream& ofs, std::vector<tap_t>& taps) {
    ofs << "$version Generated by Vortex Scope Analyzer $end" << std::endl;
    ofs << "$timescale 1 ns $end" << std::endl; 
    ofs << "$scope module TOP $end" << std::endl;
    ofs << " $var reg 1 0 clk $end" << std::endl;

    std::unordered_map<std::string, std::unordered_set<std::string>> hierarchy;
    std::unordered_set<std::string> heads;
    std::unordered_map<std::string, tap_t*> tails;

    // Build hierarchy
    for (auto& tap : taps) {
        std::vector<std::string> tokens = split(tap.path, '.');
        for (size_t i = 1; i < tokens.size(); ++i) {
            hierarchy[tokens[i-1]].insert(tokens[i]);
        }
        auto h = tokens[0];
        auto t = tokens[tokens.size()-1];
        heads.insert(h);
        tails[t] = &tap;
    }

    // Dump module huierarchy
    for (auto& head : heads) {
        dump_module(ofs, head, hierarchy, tails, 1);
    }

    ofs << "$upscope $end" << std::endl;    
    ofs << "enddefinitions $end" << std::endl;
}

static tap_t* find_nearest_tap(std::vector<tap_t>& taps) {
    tap_t* nearest = nullptr;
    for (auto& tap : taps) {
        if (tap.cur_frame == tap.frames)
            continue;
        if (nearest != nullptr) {
            if (tap.ticks < nearest->ticks)
                nearest = &tap;                
        } else {
            nearest = &tap;
        }
    }
    return nearest;
}

static int dump_tap(std::ofstream& ofs, tap_t* tap, vx_device_h hdevice) {
    uint32_t signal_offset = 0;   
    uint32_t frame_offset = 0;
    uint64_t word;

    std::vector<char> signal_data(tap->width);
    auto signal_it = tap->signals.rbegin();
    uint32_t signal_width = signal_it->width;

    do {
        // read data
        uint64_t cmd_data = (tap->id << 3) | CMD_GET_DATA;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_data));        
        CHECK_ERR(g_callback.registerRead(hdevice, &word));        
        do {            
            uint32_t word_offset = frame_offset % 64;
            signal_data[signal_width - signal_offset - 1] = ((word >> word_offset) & 0x1) ? '1' : '0';
            ++signal_offset;
            ++frame_offset;
            if (signal_offset == signal_width) {
                signal_data[signal_width] = 0; // string null termination
                ofs << 'b' << signal_data.data() << ' ' << signal_it->id << std::endl;
                if (frame_offset == tap->width) {
                    // end-of-frame
                    ++tap->cur_frame;
                    if (tap->cur_frame != tap->frames) {
                        // read next delta
                        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_data));      
                        CHECK_ERR(g_callback.registerRead(hdevice, &word));
                        tap->ticks += word;                        
                
                        if (0 == (tap->cur_frame % FRAME_FLUSH_SIZE)) {
                            ofs << std::flush;
                            std::cout << std::dec << "[SCOPE] flush tap #" << tap->id << ": "<< tap->cur_frame << "/" << tap->frames << " frames" << std::endl;
                        }
                    }
                    break; 
                }
                signal_offset = 0;
                ++signal_it;
                signal_width = signal_it->width;
            }
        } while ((frame_offset % 64) != 0);
    } while (frame_offset != tap->width);
    return 0;
}

int vx_scope_start(scope_callback_t* callback, vx_device_h hdevice, uint64_t start_time, uint64_t stop_time) {    
    if (nullptr == hdevice || nullptr == callback)
        return -1;

    const char* json_path = getenv("SCOPE_JSON_PATH");
    std::ifstream ifs(json_path);
    if (!ifs) {
        std::cerr << "[SCOPE] error: cannot open scope manifest file: " << json_path << std::endl;
        return -1;
    }
    auto json_obj = json::parse(ifs);
    if (json_obj.is_null()) {
        std::cerr << "[SCOPE] error: invalid scope manifest file: " << json_path << std::endl;
        return -1;
    }

    g_callback = *callback;   

    // validate scope manifest
    for (auto& tap : json_obj["taps"]) {
        auto id = tap["id"].get<uint32_t>();
        auto width = tap["width"].get<uint32_t>();
        
        uint64_t cmd_width = (id << 3) | CMD_GET_WIDTH;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_width));
        uint64_t dev_width;
        CHECK_ERR(g_callback.registerRead(hdevice, &dev_width));
        if (width != dev_width) {
            std::cerr << "[SCOPE] error: invalid tap #" << id << " width, actual=" << dev_width << ", expected=" << width << std::endl;
            return 1;
        }
    }

    // set stop time
    if (stop_time != uint64_t(-1)) {
        std::cout << "[SCOPE] stop time: " << std::dec << stop_time << "s" << std::endl;
        for (auto& tap : json_obj["taps"]) {
            auto id = tap["id"].get<uint32_t>();
            uint64_t cmd_stop = (stop_time << 11) | (id << 3) | CMD_SET_STOP;
            CHECK_ERR(g_callback.registerWrite(hdevice, cmd_stop));
        }        
    }

    // start recording
    if (start_time != uint64_t(-1)) {  
        std::cout << "[SCOPE] start time: " << std::dec << start_time << "s" << std::endl;
        for (auto& tap : json_obj["taps"]) {
            auto id = tap["id"].get<uint32_t>();
            uint64_t cmd_start = (start_time << 11) | (id << 3) | CMD_SET_START;
            CHECK_ERR(g_callback.registerWrite(hdevice, cmd_start));
        }        
    }

    return 0;
}

int vx_scope_stop(vx_device_h hdevice) {
    if (nullptr == hdevice)
        return -1;

    std::vector<tap_t> taps;

    {
        const char* json_path = getenv("SCOPE_JSON_PATH");
        std::ifstream ifs(json_path);
        auto json_obj = json::parse(ifs);
        if (json_obj.is_null())
            return 0;

        uint32_t signal_id = 1;

        for (auto& tap : json_obj["taps"]) {
            tap_t _tap;
            _tap.id    = tap["id"].get<uint32_t>();
            _tap.width = tap["width"].get<uint32_t>();
            _tap.path  = tap["path"].get<std::string>();
            _tap.ticks = 0;
            _tap.frames = 0;
            _tap.cur_frame = 0;            

            for (auto& signal : tap["signals"]) {
                auto name  = signal[0].get<std::string>();
                auto width = signal[1].get<uint32_t>();
                _tap.signals.push_back({signal_id, name, width});
                ++signal_id;
            }

            taps.emplace_back(std::move(_tap));
        }
    }

    // stop recording
    for (auto& tap : taps) {
        uint64_t cmd_stop = (0 << 11) | (tap.id << 3) | CMD_SET_STOP;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_stop));
    }

    std::cout << "[SCOPE] trace dump begin..." << std::endl;

    std::ofstream ofs("scope.vcd");

    dump_header(ofs, taps);

    // load trace info
    for (auto& tap : taps) {
        uint64_t count, start, delta;

        // get count
        uint64_t cmd_count = (tap.id << 3) | CMD_GET_COUNT;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_count));
        CHECK_ERR(g_callback.registerRead(hdevice, &count));   

        // get start    
        uint64_t cmd_start = (tap.id << 3) | CMD_GET_START;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_start));
        CHECK_ERR(g_callback.registerRead(hdevice, &start));

        // get data
        uint64_t cmd_data = (tap.id << 3) | CMD_GET_DATA;
        CHECK_ERR(g_callback.registerWrite(hdevice, cmd_data));
        CHECK_ERR(g_callback.registerRead(hdevice, &delta));

        tap.frames = count;
        tap.ticks = start + delta;

        std::cout << std::dec << "[SCOPE] tap #" << tap.id << ": width=" << tap.width << ", num_frames=" << tap.frames << ", start_time=" << tap.ticks << ", path=" << tap.path << std::endl;
    }  

    uint64_t timestamp = 0;

    while (true) {
        // find the nearest tap
        auto tap = find_nearest_tap(taps);
        if (tap == nullptr)
            break;
        // advance clock
        timestamp = dump_clock(ofs, tap->ticks + 1, timestamp);        
        // dump tap
        CHECK_ERR(dump_tap(ofs, tap, hdevice));
    };

    std::cout << "[SCOPE] trace dump done! - " << (timestamp/2) << " cycles" << std::endl;

    return 0;
}
