#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <list>

#include "sim.h"

using namespace std;

// Helper function for robust integer log2 to avoid floating point issues
uint32_t int_log2(uint32_t n) {
    if (n == 0) return 0;
    uint32_t count = 0;
    while (n > 1) {
        n >>= 1;
        count++;
    }
    return count;
}

//////////////////////////////////////////////////
// Cache Class Method Implementations
//////////////////////////////////////////////////

Cache::Cache(uint32_t size, uint32_t assoc, uint32_t blocksize, std::string cache_name, uint32_t pref_n, uint32_t pref_m) {
    this->name = cache_name;
    this->assoc = assoc;
    this->blocksize = blocksize;
    this->pref_n = pref_n;
    this->pref_m = pref_m;

    uint32_t num_blocks = size / blocksize;
    this->sets = num_blocks / assoc;

    this->block_offset_bits = int_log2(blocksize);
    this->index_bits = int_log2(this->sets);
    this->tag_bits = 32 - this->block_offset_bits - this->index_bits;

    cache_data.resize(this->sets);
    for (uint32_t i = 0; i < this->sets; ++i) {
        cache_data[i].resize(this->assoc);
        for (uint32_t j = 0; j < this->assoc; ++j) {
            cache_data[i][j].lru_counter = j;
        }
    }

    if (this->pref_n > 0) {
        stream_buffers.resize(this->pref_n);
        for (uint32_t i = 0; i < this->pref_n; ++i) {
            stream_buffers[i].block_valid.resize(this->pref_m, false);
            sb_lru_list.push_back(i);
        }
        last_miss_block_addr = 0xFFFFFFFF; // Initialize to invalid
    }
}

uint32_t Cache::get_tag(uint32_t address) {
    return address >> (block_offset_bits + index_bits);
}

uint32_t Cache::get_index(uint32_t address) {
    if (index_bits == 0) return 0;
    uint32_t mask = (1 << index_bits) - 1;
    return (address >> block_offset_bits) & mask;
}

void Cache::update_lru(uint32_t set_index, int mru_way) {
    uint32_t old_lru = cache_data[set_index][mru_way].lru_counter;
    if (old_lru == 0) return;

    for (uint32_t i = 0; i < assoc; ++i) {
        if (cache_data[set_index][i].lru_counter < old_lru) {
            cache_data[set_index][i].lru_counter++;
        }
    }
    cache_data[set_index][mru_way].lru_counter = 0;
}

int Cache::find_lru_victim(uint32_t set_index) {
    for (uint32_t i = 0; i < assoc; ++i) {
        if (cache_data[set_index][i].lru_counter == (assoc - 1)) {
            return i;
        }
    }
    return 0; // Should not be reached in a full set
}

// Helper function to check if a block is in a stream buffer
bool Cache::is_block_in_stream_buffer(uint32_t block_addr, int sb_idx) {
    if (!stream_buffers[sb_idx].valid) return false;
    if (block_addr < stream_buffers[sb_idx].start_block_addr || 
        block_addr > stream_buffers[sb_idx].end_block_addr) return false;
    
    uint32_t offset = block_addr - stream_buffers[sb_idx].start_block_addr;
    return offset < stream_buffers[sb_idx].block_valid.size() && 
           stream_buffers[sb_idx].block_valid[offset];
}

// Helper function to allocate a block in cache (used by prefetcher)
void Cache::allocate_block_in_cache(uint32_t address) {
    uint32_t tag = get_tag(address);
    uint32_t index = get_index(address);
    
    // Check if already in cache
    for (uint32_t way = 0; way < assoc; ++way) {
        if (cache_data[index][way].valid && cache_data[index][way].tag == tag) {
            return; // Already in cache
        }
    }
    
    // Find victim and allocate
    int victim_way = -1;
    for (uint32_t j = 0; j < assoc; ++j) {
        if (!cache_data[index][j].valid) {
            victim_way = j;
            break;
        }
    }
    if (victim_way == -1) {
        victim_way = find_lru_victim(index);
    }
    
    if (cache_data[index][victim_way].valid && cache_data[index][victim_way].dirty) {
        writebacks++;
        uint32_t victim_tag = cache_data[index][victim_way].tag;
        uint32_t victim_address = (victim_tag << (index_bits + block_offset_bits)) | (index << block_offset_bits);
        if (next_level != nullptr) {
            next_level->access(victim_address, 'w');
        }
    }
    
    cache_data[index][victim_way].valid = true;
    cache_data[index][victim_way].tag = tag;
    cache_data[index][victim_way].dirty = false;
    update_lru(index, victim_way);
}

void Cache::print_stream_buffers() {
    if (pref_n == 0) return;
    
    printf("===== Stream Buffer(s) contents =====\n");
    for (uint32_t i = 0; i < pref_n; i++) {
        if (stream_buffers[i].valid) {
            for (uint32_t j = 0; j < pref_m; j++) {
                if (stream_buffers[i].block_valid[j]) {
                    uint32_t block_addr = stream_buffers[i].start_block_addr + j;
                    printf(" %x", block_addr);
                }
            }
            printf(" \n");
        }
    }
}

void Cache::handle_prefetch(uint32_t miss_address) {
    if (pref_n == 0) return;

    uint32_t miss_block_addr = miss_address / blocksize;
    
    // Scenario #1: Miss in both cache and stream buffers - create new stream
    int hit_sb_index = sb_lru_list.back(); // Get LRU
    sb_lru_list.pop_back();
    sb_lru_list.push_front(hit_sb_index); // Move to MRU
    
    // Initialize new stream buffer
    stream_buffers[hit_sb_index].valid = true;
    stream_buffers[hit_sb_index].start_block_addr = miss_block_addr + 1;
    stream_buffers[hit_sb_index].end_block_addr = miss_block_addr + pref_m;
    
    // Clear all blocks initially
    for (uint32_t i = 0; i < pref_m; i++) {
        stream_buffers[hit_sb_index].block_valid[i] = false;
    }
    
    // Prefetch M blocks into the stream buffer
    for (uint32_t i = 0; i < pref_m; i++) {
        uint32_t prefetch_block_addr = miss_block_addr + 1 + i;
        uint32_t prefetch_addr = prefetch_block_addr * blocksize;
        prefetches++;
        
        // Fetch from next level (memory)
        if (next_level != nullptr) {
            next_level->access(prefetch_addr, 'r');
        }
        
        // Mark as valid in stream buffer
        stream_buffers[hit_sb_index].block_valid[i] = true;
    }
}

void Cache::access(uint32_t address, char rw, bool is_prefetch) {
    if (is_prefetch) {
        reads_prefetch++;
    } else {
        if (rw == 'r') reads++;
        else writes++;
    }

    uint32_t tag = get_tag(address);
    uint32_t index = get_index(address);
    uint32_t block_addr = address / blocksize;
    int hit_way = -1;

    // Check cache first
    for (uint32_t way = 0; way < assoc; ++way) {
        if (cache_data[index][way].valid && cache_data[index][way].tag == tag) {
            hit_way = way;
            break;
        }
    }

    if (hit_way != -1) { // Hit in cache
        update_lru(index, hit_way);
        if (rw == 'w' && !is_prefetch) {
            cache_data[index][hit_way].dirty = true;
        }
        
        // Check stream buffers for hits (scenarios 3 and 4)
        if (!is_prefetch) {
            // Search in MRU order (front to back) to find the most-recently-used hit
            for (auto it = sb_lru_list.begin(); it != sb_lru_list.end(); ++it) {
                int sb_idx = *it;
                if (is_block_in_stream_buffer(block_addr, sb_idx)) {
                    // Hit in stream buffer - continue stream (Scenario #4)
                    sb_lru_list.erase(it);
                    sb_lru_list.push_front(sb_idx); // Move to MRU
                    
                    // Find the offset of the hit block
                    uint32_t offset = block_addr - stream_buffers[sb_idx].start_block_addr;
                    
                    // Remove all blocks before and including the hit block
                    for (uint32_t i = 0; i <= offset; i++) {
                        stream_buffers[sb_idx].block_valid[i] = false;
                    }
                    
                    // Shift remaining blocks up
                    for (uint32_t i = 0; i < pref_m - offset - 1; i++) {
                        stream_buffers[sb_idx].block_valid[i] = stream_buffers[sb_idx].block_valid[i + offset + 1];
                    }
                    
                    // Clear the shifted positions
                    for (uint32_t i = pref_m - offset - 1; i < pref_m; i++) {
                        stream_buffers[sb_idx].block_valid[i] = false;
                    }
                    
                    // Update stream buffer range - shift by the number of consumed blocks
                    stream_buffers[sb_idx].start_block_addr += offset + 1;
                    stream_buffers[sb_idx].end_block_addr += offset + 1;
                    
                    // Prefetch new blocks to fill the empty slots
                    for (uint32_t i = 0; i <= offset; i++) {
                        uint32_t new_prefetch_block_addr = stream_buffers[sb_idx].end_block_addr - offset + i;
                        uint32_t new_prefetch_addr = new_prefetch_block_addr * blocksize;
                        prefetches++;
                        
                        // Fetch from next level (memory)
                        if (next_level != nullptr) {
                            next_level->access(new_prefetch_addr, 'r');
                        }
                        
                        // Mark as valid in stream buffer
                        stream_buffers[sb_idx].block_valid[pref_m - offset - 1 + i] = true;
                    }
                    break;
                }
            }
        }
    } else { // Miss in cache
        // Check stream buffers for hits (scenarios 1 and 2)
        bool hit_in_stream_buffer = false;
        if (!is_prefetch) {
            // Search in MRU order (front to back) to find the most-recently-used hit
            for (auto it = sb_lru_list.begin(); it != sb_lru_list.end(); ++it) {
                int sb_idx = *it;
                if (is_block_in_stream_buffer(block_addr, sb_idx)) {
                    // Hit in stream buffer - copy to cache and continue stream (Scenario #2)
                    sb_lru_list.erase(it);
                    sb_lru_list.push_front(sb_idx); // Move to MRU
                    
                    // Copy block from stream buffer to cache
                    allocate_block_in_cache(address);
                    
                    // Find the offset of the hit block
                    uint32_t offset = block_addr - stream_buffers[sb_idx].start_block_addr;
                    
                    // Remove all blocks before and including the hit block
                    for (uint32_t i = 0; i <= offset; i++) {
                        stream_buffers[sb_idx].block_valid[i] = false;
                    }
                    
                    // Shift remaining blocks up
                    for (uint32_t i = 0; i < pref_m - offset - 1; i++) {
                        stream_buffers[sb_idx].block_valid[i] = stream_buffers[sb_idx].block_valid[i + offset + 1];
                    }
                    
                    // Clear the shifted positions
                    for (uint32_t i = pref_m - offset - 1; i < pref_m; i++) {
                        stream_buffers[sb_idx].block_valid[i] = false;
                    }
                    
                    // Update stream buffer range - shift by the number of consumed blocks
                    stream_buffers[sb_idx].start_block_addr += offset + 1;
                    stream_buffers[sb_idx].end_block_addr += offset + 1;
                    
                    // Prefetch new blocks to fill the empty slots
                    for (uint32_t i = 0; i <= offset; i++) {
                        uint32_t new_prefetch_block_addr = stream_buffers[sb_idx].end_block_addr - offset + i;
                        uint32_t new_prefetch_addr = new_prefetch_block_addr * blocksize;
                        prefetches++;
                        
                        // Fetch from next level (memory)
                        if (next_level != nullptr) {
                            next_level->access(new_prefetch_addr, 'r');
                        }
                        
                        // Mark as valid in stream buffer
                        stream_buffers[sb_idx].block_valid[pref_m - offset - 1 + i] = true;
                    }
                    hit_in_stream_buffer = true;
                    break;
                }
            }
        }
        
        if (!hit_in_stream_buffer) {
            // Miss in both cache and stream buffers
            if (is_prefetch) {
                read_misses_prefetch++;
            } else {
                if (rw == 'r') read_misses++;
                else write_misses++;
            }

            // Trigger prefetcher on demand misses (only when miss in both cache and stream buffers)
            if (!is_prefetch) {
                handle_prefetch(address);
            }

            if (next_level != nullptr) {
                next_level->access(address, 'r');
            }

            int victim_way = -1;
            for (uint32_t i = 0; i < assoc; ++i) {
                if (!cache_data[index][i].valid) {
                    victim_way = i;
                    break;
                }
            }
            if (victim_way == -1) {
                victim_way = find_lru_victim(index);
            }

            if (cache_data[index][victim_way].valid && cache_data[index][victim_way].dirty) {
                writebacks++;
                uint32_t victim_tag = cache_data[index][victim_way].tag;
                uint32_t victim_address = (victim_tag << (index_bits + block_offset_bits)) | (index << block_offset_bits);
                if (next_level != nullptr) {
                    next_level->access(victim_address, 'w');
                }
            }

            cache_data[index][victim_way].valid = true;
            cache_data[index][victim_way].tag = tag;
            cache_data[index][victim_way].dirty = (rw == 'w' && !is_prefetch);
            update_lru(index, victim_way);
        }
    }
}

void Cache::print_contents() {
    cout << "===== " << name << " contents =====" << endl;
    for (uint32_t i = 0; i < sets; ++i) {
        cout << "set " << dec << i << ": ";
        for (uint32_t j = 0; j < assoc; ++j) {
            if (cache_data[i][j].valid) {
                cout << hex << cache_data[i][j].tag << (cache_data[i][j].dirty ? " D" : "  ") << " ";
            }
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 9) {
        printf("Error: Expected 8 command-line arguments but was provided %d.\n", (argc - 1));
        exit(EXIT_FAILURE);
    }
    
    cache_params_t params;
    params.BLOCKSIZE = (uint32_t)atoi(argv[1]);
    params.L1_SIZE = (uint32_t)atoi(argv[2]);
    params.L1_ASSOC = (uint32_t)atoi(argv[3]);
    params.L2_SIZE = (uint32_t)atoi(argv[4]);
    params.L2_ASSOC = (uint32_t)atoi(argv[5]);
    params.PREF_N = (uint32_t)atoi(argv[6]);
    params.PREF_M = (uint32_t)atoi(argv[7]);
    char* trace_file = argv[8];

    printf("===== Simulator configuration =====\n");
    printf("BLOCKSIZE: %u\n", params.BLOCKSIZE);
    printf("L1_SIZE: %u\n", params.L1_SIZE);
    printf("L1_ASSOC: %u\n", params.L1_ASSOC);
    printf("L2_SIZE: %u\n", params.L2_SIZE);
    printf("L2_ASSOC: %u\n", params.L2_ASSOC);
    printf("PREF_N: %u\n", params.PREF_N);
    printf("PREF_M: %u\n", params.PREF_M);
    printf("trace_file: %s\n", trace_file);
    printf("\n");

    Cache* l1_cache = nullptr;
    Cache* l2_cache = nullptr;

    if (params.L1_SIZE > 0) {
        // Prefetcher goes on the last level cache (L1 if no L2, L2 if L2 exists)
        uint32_t l1_pref_n = (params.L2_SIZE == 0) ? params.PREF_N : 0;
        uint32_t l1_pref_m = (params.L2_SIZE == 0) ? params.PREF_M : 0;
        l1_cache = new Cache(params.L1_SIZE, params.L1_ASSOC, params.BLOCKSIZE, "L1", l1_pref_n, l1_pref_m);
    }
    if (params.L2_SIZE > 0) {
        l2_cache = new Cache(params.L2_SIZE, params.L2_ASSOC, params.BLOCKSIZE, "L2", params.PREF_N, params.PREF_M);
    }
    if (l1_cache && l2_cache) {
        l1_cache->next_level = l2_cache;
    }

    FILE* fp = fopen(trace_file, "r");
    if (fp == NULL) {
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }

    char rw;
    uint32_t addr;
    while (fscanf(fp, " %c %x", &rw, &addr) == 2) {
        if (l1_cache) {
            l1_cache->access(addr, rw);
        }
    }
    fclose(fp);

    if (l1_cache) l1_cache->print_contents();
    if (l2_cache) l2_cache->print_contents();
    
    // Print stream buffer contents if any cache has prefetcher
    if (l1_cache && l1_cache->pref_n > 0) {
        l1_cache->print_stream_buffers();
    } else if (l2_cache && l2_cache->pref_n > 0) {
        l2_cache->print_stream_buffers();
    }
    
    cout << endl;

    uint32_t l1_reads = l1_cache ? l1_cache->reads : 0;
    uint32_t l1_read_misses = l1_cache ? l1_cache->read_misses : 0;
    uint32_t l1_writes = l1_cache ? l1_cache->writes : 0;
    uint32_t l1_write_misses = l1_cache ? l1_cache->write_misses : 0;
    uint32_t l1_writebacks = l1_cache ? l1_cache->writebacks : 0;
    uint32_t l1_prefetches = l1_cache ? l1_cache->prefetches : 0;

    uint32_t l2_reads_demand = l2_cache ? l2_cache->reads : 0;
    uint32_t l2_read_misses_demand = l2_cache ? l2_cache->read_misses : 0;
    uint32_t l2_reads_prefetch = l2_cache ? l2_cache->reads_prefetch : 0;
    uint32_t l2_read_misses_prefetch = l2_cache ? l2_cache->read_misses_prefetch : 0;
    uint32_t l2_writes = l2_cache ? l2_cache->writes : 0;
    uint32_t l2_write_misses = l2_cache ? l2_cache->write_misses : 0;
    uint32_t l2_writebacks = l2_cache ? l2_cache->writebacks : 0;
    uint32_t l2_prefetches = l2_cache ? l2_cache->prefetches : 0;

    uint32_t mem_traffic = 0;
    if (l2_cache) {
        mem_traffic = l2_read_misses_demand + l2_read_misses_prefetch + l2_write_misses + l2_writebacks + l2_prefetches;
    } else if (l1_cache) {
        mem_traffic = l1_read_misses + l1_write_misses + l1_writebacks + l1_prefetches;
    }

    printf("===== Measurements =====\n");
    printf("a. L1 reads:                   %u\n", l1_reads);
    printf("b. L1 read misses:             %u\n", l1_read_misses);
    printf("c. L1 writes:                  %u\n", l1_writes);
    printf("d. L1 write misses:            %u\n", l1_write_misses);
    double l1_miss_rate = (l1_reads + l1_writes == 0) ? 0 : (double)(l1_read_misses + l1_write_misses) / (l1_reads + l1_writes);
    printf("e. L1 miss rate:               %.4f\n", l1_miss_rate);
    printf("f. L1 writebacks:              %u\n", l1_writebacks);
    printf("g. L1 prefetches:              %u\n", l1_prefetches);
    printf("h. L2 reads (demand):          %u\n", l2_reads_demand);
    printf("i. L2 read misses (demand):    %u\n", l2_read_misses_demand);
    printf("j. L2 reads (prefetch):        %u\n", l2_reads_prefetch);
    printf("k. L2 read misses (prefetch):  %u\n", l2_read_misses_prefetch);
    printf("l. L2 writes:                  %u\n", l2_writes);
    printf("m. L2 write misses:            %u\n", l2_write_misses);
    double l2_miss_rate = (l2_reads_demand == 0) ? 0 : (double)l2_read_misses_demand / l2_reads_demand;
    printf("n. L2 miss rate:               %.4f\n", l2_miss_rate);
    printf("o. L2 writebacks:              %u\n", l2_writebacks);
    printf("p. L2 prefetches:              %u\n", l2_prefetches);
    printf("q. memory traffic:             %u\n", mem_traffic);

    delete l1_cache;
    if (l2_cache) delete l2_cache;
    
    return 0;
}