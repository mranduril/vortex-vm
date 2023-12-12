#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <future>
#include <chrono>
#include <bitset>
#include <unistd.h>

#include <vortex.h>
#include <utils.h>
#include <malloc.h>

#include <VX_config.h>
#include <VX_types.h>

#include <util.h>

#include <processor.h>
#include <arch.h>
#include <mem.h>
#include <constants.h>

#ifndef NDEBUG
#define DBGPRINT(format, ...) do { printf("[VXDRV] " format "", ##__VA_ARGS__); } while (0)
#else
#define DBGPRINT(format, ...) ((void)0)
#endif

using namespace vortex;

uint64_t bits(uint64_t addr, uint8_t s_idx, uint8_t e_idx)
{
    return (addr >> s_idx) & ((1 << (e_idx - s_idx + 1)) - 1);
}
bool bit(uint64_t addr, uint8_t idx)
{
    return (addr) & (1 << idx);
}

///////////////////////////////////////////////////////////////////////////////

class vx_device;

class vx_buffer {
public:
    vx_buffer(uint64_t size, vx_device* device) 
        : size_(size)
        , device_(device) {
        uint64_t aligned_asize = aligned_size(size, CACHE_BLOCK_SIZE);
        data_ = aligned_malloc(aligned_asize, CACHE_BLOCK_SIZE);
        // set uninitialized data to "baadf00d"
        for (uint32_t i = 0; i < aligned_asize; ++i) {
            ((uint8_t*)data_)[i] = (0xbaadf00d >> ((i & 0x3) * 8)) & 0xff;
        }
    }

    ~vx_buffer() {
        if (data_) {
            aligned_free(data_);
        }
    }

    void* data() const {
        return data_;
    }

    uint64_t size() const {
        return size_;
    }

    vx_device* device() const {
        return device_;
    }

private:
    uint64_t   size_;
    vx_device* device_;
    void*      data_;
};

///////////////////////////////////////////////////////////////////////////////

class vx_device {    
public:
    vx_device() 
        : arch_(NUM_CORES * NUM_CLUSTERS, NUM_WARPS, NUM_THREADS)
        , ram_(RAM_PAGE_SIZE)
        , processor_(arch_)
        , mem_allocator_(
            ALLOC_BASE_ADDR,
            ALLOC_MAX_ADDR,
            RAM_PAGE_SIZE,
            CACHE_BLOCK_SIZE) 
    {
        // attach memory module
        processor_.attach_ram(&ram_);

        //Sets more
        set_processor_satp(VM_ADDR_MODE);
    }

    ~vx_device() {
        if (future_.valid()) {
            future_.wait();
        }
    }    

    int map_local_mem(uint64_t size, uint64_t* dev_maddr) 
    {
        bool skip = false;
        if (*dev_maddr == STARTUP_ADDR || *dev_maddr == 0x7FFFF000) {
            skip = true;
        }

        if (get_mode() == VA_MODE::BARE) {
            return 0;
        }

        uint64_t ppn = *dev_maddr >> 12;
        uint64_t init_pAddr = *dev_maddr;
        uint64_t init_vAddr = *dev_maddr + 0xf0000000; // vpn will change, but we want to return the vpn of the beginning of the virtual allocation
        init_vAddr = (init_vAddr >> 12) << 12;
        uint64_t vpn;

        //dev_maddr can be of size greater than a page, but we have to map and update
        //page tables on a page table granularity. So divide the allocation into pages.
        for (ppn = (*dev_maddr) >> 12; ppn < ((*dev_maddr) >> 12) + (size/RAM_PAGE_SIZE) + 1; ppn++)
        {
            //Currently a 1-1 mapping is used, this can be changed here to support different
            //mapping schemes
            vpn = skip ? ppn : ppn + 0xf0000;
            //vpn = ppn;

            //If ppn to vpn mapping doesnt exist.
            if (addr_mapping.find(vpn) == addr_mapping.end())
            {
                //Create mapping.
                update_page_table(ppn, vpn);
                addr_mapping[vpn] = ppn;
            }
        }

        uint64_t size_bits;
        if (skip) {
            return 0;
        }

        *dev_maddr = init_vAddr; // commit vpn to be returned to host
        
        return 0;
    }

    int mem_alloc(uint64_t size, uint64_t* dev_maddr) {
        int err = mem_allocator_.allocate(size, dev_maddr);
        map_local_mem(size, dev_maddr);
        return err;
    }

    int mem_free(uint64_t dev_maddr) {
        return mem_allocator_.release(dev_maddr);
    }

    uint64_t mem_used() const {
        return mem_allocator_.allocated();
    }

    int upload(const void* src, uint64_t dest_addr, uint64_t size, uint64_t src_offset) {
        uint64_t asize = aligned_size(size, CACHE_BLOCK_SIZE);
        uint64_t pAddr = dest_addr; // map_local_mem overwrites the provided dest_addr, so store away physical destination address
        if (dest_addr + asize > LOCAL_MEM_SIZE)
            return -1;

        if (dest_addr >= STARTUP_ADDR) {
            map_local_mem(asize,&dest_addr);
        } else if (dest_addr >= 0x7fff0000)
        {
            map_local_mem(asize,&dest_addr);
        }
        ram_.write((const uint8_t*)src + src_offset, pAddr, asize);        
        return 0;
    }

    int download(void* dest, uint64_t src_addr, uint64_t size, uint64_t dest_offset) {
        uint64_t asize = aligned_size(size, CACHE_BLOCK_SIZE);
        if (src_addr + asize > LOCAL_MEM_SIZE)
            return -1;

        ram_.read((uint8_t*)dest + dest_offset, src_addr, asize);
        
        /*
        printf("VXDRV: download %d bytes from 0x%x\n", size, src_addr);
        for (int i = 0; i < size; i += 4) {
            printf("mem-read: 0x%x -> 0x%x\n", src_addr + i, *(uint32_t*)((uint8_t*)dest + dest_offset + i));
        }*/
        
        return 0;
    }

    int start() {  
        // ensure prior run completed
        if (future_.valid()) {
            future_.wait();
        }
        
        // start new run
        future_ = std::async(std::launch::async, [&]{
            processor_.run();
        });
        
        return 0;
    }

    int wait(uint64_t timeout) {
        if (!future_.valid())
            return 0;
        uint64_t timeout_sec = timeout / 1000;
        std::chrono::seconds wait_time(1);
        for (;;) {
            // wait for 1 sec and check status
            auto status = future_.wait_for(wait_time);
            if (status == std::future_status::ready 
             || 0 == timeout_sec--)
                break;
        }
        return 0;
    } 

    void set_processor_satp(VA_MODE mode)
    {
        uint32_t satp;
        if (mode == VA_MODE::BARE)
            satp = 0;
        else if (mode == VA_MODE::SV64)
        {
            satp = alloc_page_table();
        }
        processor_.set_satp(satp);
    }

    uint32_t get_ptbr()
    {        

        return processor_.get_satp();
    }

    VA_MODE get_mode()
    {
        return VA_MODE::SV64;
    }  

    void update_page_table(uint64_t pAddr, uint64_t vAddr) {
        std::cout << "mapping PPN 0x" << std::hex << pAddr << " to VPN 0x" << std::hex << vAddr << ":" << std::endl;
        std::cout << "\t";
        vAddr = vAddr << 12;
        //Updating page table with the following mapping of (vAddr) to (pAddr).
        uint64_t ppn_1, pte_addr, pte_bytes;
        uint64_t vpn_1 = bits(vAddr, 22, 31);
        uint64_t vpn_0 = bits(vAddr, 12, 21);

        //Read first level PTE.
        pte_addr = get_ptbr() + vpn_1 * PTE_SIZE;
        pte_bytes = read_pte(pte_addr);     

        if ( bit(pte_bytes, 0) && ((pte_bytes & 0xFFFFFFFF) != 0xbaadf00d))
        {
            //If valid bit set, proceed to next level using new ppn form PTE.
            ppn_1 = (pte_bytes >> 10);
            std::cout << pte_bytes;
        }
        else
        {
            //If valid bit not set, allocate a second level page table
            // in device memory and store ppn in PTE. Set rwx = 000 in PTE
            //to indicate this is a pointer to the next level of the page table.
            ppn_1 = alloc_page_table();
            pte_bytes = ( (ppn_1 << 10) | 0b0000000001);
            std::cout << pte_bytes;
            write_pte(pte_addr, pte_bytes);
        }
        std::cout << " --> " << std::endl << "\t\t";

        //Read second level PTE.
        pte_addr = ppn_1 + vpn_0 * PTE_SIZE;
        pte_bytes = read_pte(pte_addr); 
    
        if ( bit(pte_bytes, 0) && ((pte_bytes & 0xFFFFFFFF) != 0xbaadf00d))
        {
            std::cout << "ERROR, shouldn't be here" << std::endl;
            //If valid bit is set, then the page is already allocated.
            //Should not reach this point, a sanity check.
        }
        else
        {
            //If valid bit not set, write ppn of pAddr in PTE. Set rwx = 111 in PTE
            //to indicate this is a leaf PTE and has the stated permissions.
            pte_bytes = ( (pAddr << 10) | 0b0000001111) ;
            std::cout << pte_bytes;
            write_pte(pte_addr, pte_bytes);

            //If super paging is enabled.
            if (SUPER_PAGING)
            {
                //Check if this second level Page Table can be promoted to a super page. Brute force 
                //method is used to iterate over all PTE entries of the table and check if they have 
                //their valid bit set.
                bool superpage = true;
                for(int i = 0; i < 1024; i++)
                {
                    pte_addr = (ppn_1 << 12) + i;
                    pte_bytes = read_pte(pte_addr); 
                  
                    if (!bit(pte_bytes, 0) && ((pte_bytes & 0xFFFFFFFF) != 0xbaadf00d))
                    {
                        superpage = false;
                        break;
                    }
                }
                if (superpage)
                {
                    //This can be promoted to a super page. Set root PTE to the first PTE of the 
                    //second level. This is because the first PTE of the second level already has the
                    //correct PPN1, PPN0 set to zero and correct access bits.
                    pte_addr = (ppn_1 << 12);
                    pte_bytes = read_pte(pte_addr);
                    pte_addr = get_ptbr() + vpn_1 * PTE_SIZE;
                    write_pte(pte_addr, pte_bytes);
                }
            }
        }
        std::cout << " --> " << std::endl << "\t\t\t0x" << std::hex << pAddr << std::endl;
    }

    std::pair<uint64_t, uint8_t> page_table_walk(uint64_t vAddr_bits, uint64_t* size_bits)
    {   
        uint64_t LEVELS = 2;
        vAddr_SV64_t vAddr(vAddr_bits);
        uint64_t pte_bytes;


        //Get base page table.
        uint64_t a = this->processor_.get_satp();
        int i = LEVELS - 1; 

        while(true)
        {

        //Read PTE.
        ram_.read(&pte_bytes, a+vAddr.vpn[i] * PTE_SIZE, sizeof(uint64_t));

        //pte_bytes &= 0x00000000FFFFFFFF;
        PTE_SV64_t pte(pte_bytes);
        
        //Check if it has invalid flag bits.
        if ( (pte.v == 0) | ( (pte.r == 0) & (pte.w == 1) ) )
        {
            throw Page_Fault_Exception("Page Fault : Attempted to access invalid entry. Entry: 0x");
        }

        if ( (pte.r == 0) & (pte.w == 0) & (pte.x == 0))
        {
            //Not a leaf node as rwx == 000
            i--;
            if (i < 0)
            {
            throw Page_Fault_Exception("Page Fault : No leaf node found.");
            }
            else
            {
            //Continue on to next level.
            a = (pte_bytes >> 10 );
            }
        }
        else
        {
            //Leaf node found, finished walking.
            a = (pte_bytes >> 10 ) << 12;
            break;
        }
        }

        PTE_SV64_t pte(pte_bytes);

        //Check RWX permissions according to access type.
        if (pte.r == 0)
        {
        throw Page_Fault_Exception("Page Fault : TYPE LOAD, Incorrect permissions.");
        }

        uint64_t pfn;
        if (i > 0)
        {
        //It is a super page.
        if (pte.ppn[0] != 0)
        {
            //Misss aligned super page.
            throw Page_Fault_Exception("Page Fault : Miss Aligned Super Page.");

        }
        else
        {
            //Valid super page.
            pfn = pte.ppn[1];
            *size_bits = 22;
        }
        }
        else
        {
        //Regular page.
        *size_bits = 12;
        pfn = a >> 12;
        }
        return std::make_pair(pfn, pte_bytes & 0xff);
    }

    uint64_t alloc_page_table() {
        uint64_t addr;
        mem_allocator_.allocate(RAM_PAGE_SIZE * 2, &addr);
        init_page_table(addr);
        return addr;
    }

    void init_page_table(uint64_t addr) {
        uint64_t asize = aligned_size(RAM_PAGE_SIZE * 2, CACHE_BLOCK_SIZE);
        uint8_t *src = new uint8_t[RAM_PAGE_SIZE * 2];
        for (uint64_t i = 0; i < RAM_PAGE_SIZE * 2; ++i) {
            src[i] = (0x00000000 >> ((i & 0x3) * 8)) & 0xff;
        }
        ram_.write((const uint8_t*)src, addr, asize);
    }

    void read_page_table(uint64_t addr) {
        uint8_t *dest = new uint8_t[RAM_PAGE_SIZE * 2];
        download(dest,  addr,  RAM_PAGE_SIZE * 2, 0);
        printf("VXDRV: download %d bytes from 0x%x\n", RAM_PAGE_SIZE * 2, addr);
        for (int i = 0; i < RAM_PAGE_SIZE * 2; i += 4) {
            printf("mem-read: 0x%x -> 0x%x\n", addr + i, *(uint64_t*)((uint8_t*)dest + i));
        }
    }

    void write_pte(uint64_t addr, uint64_t value = 0xbaadf00d) {
        uint8_t *src = new uint8_t[PTE_SIZE];
        for (uint64_t i = 0; i < PTE_SIZE; ++i) {
            //(value >> ((i & 0x3) * 8)) & 0xff;
            src[i] = (value >> (i * 8)) & 0xff;
        }
        ram_.write((const uint8_t*)src, addr, PTE_SIZE);
    }

    uint64_t read_pte(uint64_t addr) {
        uint8_t *dest = new uint8_t[PTE_SIZE];
        ram_.read((uint8_t*)dest, addr, PTE_SIZE);
        return *(uint64_t*)((uint8_t*)dest);
    }

    int write_dcr(uint32_t addr, uint32_t value) {
        if (future_.valid()) {
            future_.wait(); // ensure prior run completed
        }        
        processor_.write_dcr(addr, value);
        dcrs_.write(addr, value);
        return 0;
    }

    uint64_t read_dcr(uint32_t addr) const {
        return dcrs_.read(addr);
    }

private:
    Arch                arch_;
    RAM                 ram_;
    Processor           processor_;
    MemoryAllocator     mem_allocator_;
    DeviceConfig        dcrs_;
    std::future<void>   future_;
    std::unordered_map<uint64_t, uint64_t> addr_mapping;
};

///////////////////////////////////////////////////////////////////////////////

extern int vx_dev_open(vx_device_h* hdevice) {
    if (nullptr == hdevice)
        return  -1;

    auto device = new vx_device();
    if (device == nullptr)
        return -1;

    int err = dcr_initialize(device);
    if (err != 0) {
        delete device;
        return err;
    }

#ifdef DUMP_PERF_STATS
    perf_add_device(device);
#endif  

    *hdevice = device;

    DBGPRINT("device creation complete!\n");

    return 0;
}

extern int vx_dev_close(vx_device_h hdevice) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);

#ifdef DUMP_PERF_STATS
    perf_remove_device(hdevice);
#endif

    delete device;

    DBGPRINT("device destroyed!\n");

    return 0;
}

extern int vx_dev_caps(vx_device_h hdevice, uint32_t caps_id, uint64_t *value) {
    if (nullptr == hdevice)
        return  -1;

    vx_device *device = ((vx_device*)hdevice);

    switch (caps_id) {
    case VX_CAPS_VERSION:
        *value = IMPLEMENTATION_ID;
        break;
    case VX_CAPS_NUM_CORES:
        *value = NUM_CORES * NUM_CLUSTERS;        
        break;
    case VX_CAPS_NUM_WARPS:
        *value = NUM_WARPS;
        break;
    case VX_CAPS_NUM_THREADS:
        *value = NUM_THREADS;
        break;
    case VX_CAPS_CACHE_LINE_SIZE:
        *value = CACHE_BLOCK_SIZE;
        break;
    case VX_CAPS_LOCAL_MEM_SIZE:
        *value = LOCAL_MEM_SIZE;
        break;
    case VX_CAPS_KERNEL_BASE_ADDR:
        *value = (uint64_t(device->read_dcr(DCR_BASE_STARTUP_ADDR1)) << 32) | 
                device->read_dcr(DCR_BASE_STARTUP_ADDR0);
        break;    
    case VX_CAPS_ISA_FLAGS:
        *value = ((uint64_t(MISA_EXT))<<32) | ((log2floor(XLEN)-4) << 30) | MISA_STD;
        break;
    default:
        std::cout << "invalid caps id: " << caps_id << std::endl;
        std::abort();
        return -1;
    }

    return 0;
}

extern int vx_mem_alloc(vx_device_h hdevice, uint64_t size, uint64_t* dev_maddr) {
    if (nullptr == hdevice 
     || nullptr == dev_maddr
     || 0 == size)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    return device->mem_alloc(size, dev_maddr);
}

extern int vx_mem_free(vx_device_h hdevice, uint64_t dev_maddr) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    uint64_t size_bits;
    std::pair<uint64_t, uint8_t> ptw_access = device->page_table_walk(dev_maddr, &size_bits);
    uint64_t pfn = ptw_access.first;
    dev_maddr = (pfn << 12) + 0x40;

    if (0 == dev_maddr)
        return 0;

    return device->mem_free(dev_maddr);
}

extern int vx_mem_info(vx_device_h hdevice, uint64_t* mem_free, uint64_t* mem_total) {
    if (nullptr == hdevice)
        return -1;

    auto device = ((vx_device*)hdevice);
    if (mem_free) {
        *mem_free = (ALLOC_MAX_ADDR - ALLOC_BASE_ADDR) - device->mem_used();
    }
    if (mem_total) {
        *mem_total = (ALLOC_MAX_ADDR - ALLOC_BASE_ADDR);
    }
    return 0;
}

extern int vx_buf_alloc(vx_device_h hdevice, uint64_t size, vx_buffer_h* hbuffer) {
    if (nullptr == hdevice 
     || 0 >= size
     || nullptr == hbuffer)
        return -1;

    vx_device *device = ((vx_device*)hdevice);

    auto buffer = new vx_buffer(size, device);
    if (nullptr == buffer->data()) {
        delete buffer;
        return -1;
    }

    *hbuffer = buffer;

    return 0;
}

extern void* vx_host_ptr(vx_buffer_h hbuffer) {
    if (nullptr == hbuffer)
        return nullptr;

    vx_buffer* buffer = ((vx_buffer*)hbuffer);

    return buffer->data();
}

extern int vx_buf_free(vx_buffer_h hbuffer) {
    if (nullptr == hbuffer)
        return -1;

    vx_buffer* buffer = ((vx_buffer*)hbuffer);

    delete buffer;

    return 0;
}

extern int vx_copy_to_dev(vx_buffer_h hbuffer, uint64_t dev_maddr, uint64_t size, uint64_t src_offset) {
    if (nullptr == hbuffer 
     || 0 >= size)
        return -1;

    auto buffer = (vx_buffer*)hbuffer;

    if (size + src_offset > buffer->size())
        return -1;

    if (!(dev_maddr == STARTUP_ADDR) && !(dev_maddr == 0x7FFFF000)) {
        uint64_t size_bits;
        std::pair<uint64_t, uint8_t> ptw_access = buffer->device()->page_table_walk(dev_maddr, &size_bits);
        uint64_t pfn = ptw_access.first;
        dev_maddr = pfn << 12;
    }
    return buffer->device()->upload(buffer->data(), dev_maddr, size, src_offset);
}

extern int vx_copy_from_dev(vx_buffer_h hbuffer, uint64_t dev_maddr, uint64_t size, uint64_t dest_offset) {
     if (nullptr == hbuffer 
      || 0 >= size)
        return -1;

    auto buffer = (vx_buffer*)hbuffer;

    if (size + dest_offset > buffer->size())
        return -1;    

    uint64_t size_bits;
    std::pair<uint64_t, uint8_t> ptw_access = buffer->device()->page_table_walk(dev_maddr, &size_bits);
    uint64_t pfn = ptw_access.first;
    return buffer->device()->download(buffer->data(), (pfn << 12), size, dest_offset);
}

extern int vx_start(vx_device_h hdevice) {
    if (nullptr == hdevice)
        return -1;    
    
    DBGPRINT("START\n");

    vx_device *device = ((vx_device*)hdevice);

    return device->start();
}

extern int vx_ready_wait(vx_device_h hdevice, uint64_t timeout) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);

    return device->wait(timeout);
}

extern int vx_dcr_write(vx_device_h hdevice, uint32_t addr, uint64_t value) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);

    // Ensure ready for new command
    if (vx_ready_wait(hdevice, -1) != 0)
        return -1;

    DBGPRINT("DCR_WRITE: addr=0x%x, value=0x%lx\n", addr, value);
  
    return device->write_dcr(addr, value);
}