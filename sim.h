#ifndef SIM_CACHE_H
#define SIM_CACHE_H




// Command-line parameters
typedef struct {
    uint32_t BLOCKSIZE;
    uint32_t L1_SIZE;
    uint32_t L1_ASSOC;
    uint32_t L2_SIZE;
    uint32_t L2_ASSOC;
    uint32_t PREF_N;
    uint32_t PREF_M;
} cache_params_t;

// Represents a single block in the cache
struct CacheBlock {
    bool valid = false;
    bool dirty = false;
    uint32_t tag = 0;
    uint32_t lru_counter = 0;
};

// Represents a single stream buffer for the prefetcher
struct StreamBuffer {
    bool valid = false;
    uint32_t start_block_addr = 0;  // First block in the stream
    uint32_t end_block_addr = 0;    // Last block in the stream
    std::vector<bool> block_valid;  // Which blocks are valid in this stream
};

class Cache {
public:
    // Constructor
    Cache(uint32_t size, uint32_t assoc, uint32_t blocksize, std::string cache_name, uint32_t pref_n, uint32_t pref_m);

    // Main access function, now with a flag for prefetch accesses
    void access(uint32_t address, char rw, bool is_prefetch = false);

    // Print final cache contents
    void print_contents();

    // Pointer to the next level in the memory hierarchy
    Cache* next_level = nullptr;

    // Public counters for final summary
    uint32_t reads = 0;
    uint32_t read_misses = 0;
    uint32_t writes = 0;
    uint32_t write_misses = 0;
    uint32_t writebacks = 0;
    // Prefetching stats
    uint32_t prefetches = 0;
    uint32_t reads_prefetch = 0;
    uint32_t read_misses_prefetch = 0;

private:
    // Cache configuration
    std::string name;
    uint32_t assoc;
    uint32_t sets;
    uint32_t blocksize;
    uint32_t block_offset_bits;
    uint32_t index_bits;
    uint32_t tag_bits;

    // The actual cache storage: a 2D vector of sets and ways
    std::vector<std::vector<CacheBlock>> cache_data;

    // Prefetcher components
    uint32_t pref_m;
    std::vector<StreamBuffer> stream_buffers;
    std::list<int> sb_lru_list; // Tracks LRU order for stream buffers
    uint32_t last_miss_block_addr; // Track last demand miss block address

    // Helper functions to get address parts
    uint32_t get_tag(uint32_t address);
    uint32_t get_index(uint32_t address);

    // LRU management
    void update_lru(uint32_t set_index, int mru_way);
    int find_lru_victim(uint32_t set_index);

    // Prefetcher management
    void handle_prefetch(uint32_t miss_address);
    bool is_block_in_stream_buffer(uint32_t block_addr, int sb_idx);
    void allocate_block_in_cache(uint32_t address);
    
public:
    uint32_t pref_n;
    void print_stream_buffers();
};

#endif // SIM_CACHE_H