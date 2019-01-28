#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "settings.h"
#include "lsm_settings.h"
#include "bloomfilter.h"
#include "sha256.h"

#define CHANNEL 8
#define WAY 8
#define CHIP ((CHANNEL)*(WAY))

#define BLOCK_PER_CHIP 8
#define PAGE_PER_BLOCK 128
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP))
#define TOTAL_PAGE ((TOTAL_BLOCK)*(PAGE_PER_BLOCK))
#define PAGE_SZ (4*(K))

#define LSB (int)(log(TOTAL_BLOCK) / log(2))
#define MASK ((1 << LSB) - 1)

#define CAPACITY ((PAGE_SZ)*(TOTAL_PAGE))
#define DATA_SET (TOTAL_PAGE)

int data_sz;
uint32_t val;

typedef struct {
    uint32_t lbanum;
    int vbit;
    int level;
    int found;
    int notfound;
    double fp;
} Dat;

// Struct-Dat follows generated sequences of their own
Dat written[DATA_SET];
Dat reading[DATA_SET]; // for statistics

int read_cnt;
int write_cnt;
int true_cnt;
int false_cnt;
int found_cnt;
int notfound_cnt;

typedef struct {
    uint32_t* page;
    uint32_t* oob;
    int empty;
    int valid;
} Block;

typedef struct {
    Block* data_blk;
} Chip;

typedef struct {
    Chip** chip_arr;
} Storage;

typedef struct {
    BF*** bfchip_arr;
    uint64_t bytes_per_chip;
    uint64_t* bytes_per_blk;
    uint64_t** bytes_per_pg;
} BFManager;

Storage storage;
BFManager** global_bf; // Management in chip-level

Block** data_block;

void bloom_init() {
    double true_p=0.0, false_p=0.0;

    // Alloc memory for all bloomfilters globally
    global_bf = (BFManager**)malloc(sizeof(BFManager*) * CHIP);
    for(int c=0; c<CHIP; c++) {
        global_bf[c] = (BFManager*)malloc(sizeof(BFManager));
        global_bf[c]->bfchip_arr = (BF***)malloc(sizeof(BF**) * BLOCK_PER_CHIP);
        global_bf[c]->bytes_per_blk = (uint64_t*)malloc(sizeof(uint64_t) * BLOCK_PER_CHIP);
        global_bf[c]->bytes_per_pg = (uint64_t**)malloc(sizeof(uint64_t*) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_bf[c]->bfchip_arr[b] = (BF**)malloc(sizeof(BF*) * PAGE_PER_BLOCK);
            global_bf[c]->bytes_per_pg[b] = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
        }
    }

    // Assign different bytes per page-level in all bloomfilters separately
    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                if(p == 0) {
                    global_bf[c]->bfchip_arr[b][p] = bf_init(1, 0.0001);
                }
                else {
                    true_p = pow(0.9, (double)1/p);
                    false_p = 1 - true_p;

                    global_bf[c]->bfchip_arr[b][p] = bf_init(1, false_p);
                    global_bf[c]->bytes_per_pg[b][p] = bf_bits(1, false_p);
                    //bf_bytes += bf_bits(1, false_p);
                }

                global_bf[c]->bytes_per_blk[b] += global_bf[c]->bytes_per_pg[b][p];
            }

            //printf("Bloomfilter %lubytes in block %d of chip %d\n", bf_bytes, b, c);
            //bf_sum_bytes += bf_bytes;
            //bf_bytes = 0;
            global_bf[c]->bytes_per_chip += global_bf[c]->bytes_per_blk[b];
        }
    }

    // Alloc memory for all chips, blocks, and pages
    storage.chip_arr = (Chip**)malloc(sizeof(Chip*) * WAY);
    for(int w=0; w<WAY; w++) {
        storage.chip_arr[w] = (Chip*)malloc(sizeof(Chip) * CHANNEL);
    }

    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
            storage.chip_arr[w][ch].data_blk = (Block*)malloc(sizeof(Block) * BLOCK_PER_CHIP);
        }
    }

    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
            storage.chip_arr[w][ch].data_blk[b].page = (uint32_t*)malloc(sizeof(uint32_t) * PAGE_PER_BLOCK);
            storage.chip_arr[w][ch].data_blk[b].oob = (uint32_t*)malloc(sizeof(uint32_t) * PAGE_PER_BLOCK);
            }
        }
    }

    // Initialize values for data blocks
    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    storage.chip_arr[w][ch].data_blk[b].page[p] = \
                        storage.chip_arr[w][ch].data_blk[b].oob[p] = 99999;
                }

                storage.chip_arr[w][ch].data_blk[b].empty = 0;
                storage.chip_arr[w][ch].data_blk[b].valid = 0;
            }
        }
    }

    // Initialize values for statistics
    for(int i=0; i<DATA_SET; i++) {
        reading[i].found = reading[i].notfound = 0;
    }
}

void bloom_destroy() {
    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                bf_free(global_bf[c]->bfchip_arr[b][p]);
            }

            free(global_bf[c]->bfchip_arr[b]);
        }
        
        free(global_bf[c]->bfchip_arr);
        free(global_bf[c]);
    }

    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                free(storage.chip_arr[w][ch].data_blk[b].page);
                free(storage.chip_arr[w][ch].data_blk[b].oob);
            }

            free(storage.chip_arr[w][ch].data_blk);
        }

        free(storage.chip_arr[w]);
    }
    free(storage.chip_arr);
}

uint32_t hashing_key(uint32_t key) {
    char* string;
    Sha256Context ctx;
    SHA256_HASH hash;
    uint32_t bytes_arr[8];
    uint32_t hashkey;
    
    string = (char*)&key;

    Sha256Initialise(&ctx);
    Sha256Update(&ctx, (unsigned char*)string, sizeof(uint32_t));
    Sha256Finalise(&ctx, &hash);

    for(int i=0; i<8; i++) {
        bytes_arr[i] = ((hash.bytes[i*4] << 24) | (hash.bytes[i*4+1] << 16) | \
                (hash.bytes[i*4+2] << 8) | (hash.bytes[i*4+3]));
    }

    hashkey = bytes_arr[0];
    for(int i=1; i<8; i++) {
        hashkey ^= bytes_arr[i];
    }

    return hashkey;
}

uint32_t generate_lba(char* req, char* pttr) {
    if(!strcmp(req, "WR")) {
        if(!strcmp(pttr, "SQ")) {
            return val;
        }
        else {
            return rand() % DATA_SET;
        }
    }
    else {
        if(!strcmp(pttr, "SQ")) {
            while(written[val].vbit != 1) {
                val++;

                if(val > data_sz) {
                    val = 0;
                }
            }
            return written[val].lbanum;
        }
        else {
            int rand_idx = rand() % DATA_SET;
            return written[rand_idx].lbanum;
        }
    }
}

void bloom_write(uint32_t lba, uint32_t value, char* pttr) {
    uint32_t pbn, key, hashkey;
    int chip, way, chnl, blk, offset;
    
    while(1) {
        pbn = (lba & MASK);
        chip = pbn / BLOCK_PER_CHIP;
        way = chip / CHANNEL;
        chnl = chip % CHANNEL;
        blk = pbn % BLOCK_PER_CHIP;
        offset = storage.chip_arr[way][chnl].data_blk[blk].empty;

        if(offset >= PAGE_PER_BLOCK) {
            lba = generate_lba("WR", pttr);
        }
        else {
            break;
        }
    }
    
    key = lba + offset;
    hashkey = hashing_key(key);
    bf_set(global_bf[chip]->bfchip_arr[blk][offset], hashkey);

    storage.chip_arr[way][chnl].data_blk[blk].page[offset] = value;
    storage.chip_arr[way][chnl].data_blk[blk].oob[offset] = lba;
    storage.chip_arr[way][chnl].data_blk[blk].empty++;

    // Record working set
    written[write_cnt].lbanum = lba;
    written[write_cnt].vbit = 1;

    if(data_sz <= lba) {
        data_sz = lba;
    }
}

void bloom_read(uint32_t lba) {
    uint32_t pbn, key, hashkey;
    int chip, way, chnl, blk, offset;
//printf("LBA %u\n", lba);
    pbn = (lba & MASK);
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;
    offset = storage.chip_arr[way][chnl].data_blk[blk].empty;

    for(int idx=offset-1; idx>=0; idx--) {
        key = lba + idx;
        hashkey = hashing_key(key);

        // Bloomfilter says true
        if(bf_check(global_bf[chip]->bfchip_arr[blk][idx], hashkey) == true) {
            if(storage.chip_arr[way][chnl].data_blk[blk].oob[idx] == lba) { // Data exists
                found_cnt++;
                reading[read_cnt].lbanum = lba;
                reading[read_cnt].level = offset - idx;
                reading[read_cnt].found++;
                return;
            }
            else { // Data doesn't exist
                notfound_cnt++;
                reading[read_cnt].notfound++;
                continue;
            }
        }
        // Bloomfilter says false
        else { // Data should not exist
            false_cnt++;
        }
    }
}

void print_stats() {
    uint64_t sum = 0;

    printf("\nFOUND LEVEL OF ALL READ REQUESTS\n");
    for(int i=0; i<DATA_SET; i++) {
        printf("req_num: %d level: %d\n", i, reading[i].level);
    }

    printf("\nFALSE POSITIVE RATE OF ALL READ REQUESTS\n");
    for(int i=0; i<DATA_SET; i++) {
        reading[i].fp = (double)reading[i].notfound / (reading[i].found + reading[i].notfound);
        printf("req_num: %d false_p: %.2lf\n", i, reading[i].fp);
    }
    
    printf("\nTEST\n");
    printf("Total read requests: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("Total false num: %d\n", false_cnt);
    printf("RAF: %.2f\n", (double)(found_cnt + notfound_cnt) / read_cnt);
    
    printf("Sum of bloomfilter assigned bytes per chip: ");
    for(int c=0; c<CHIP; c++) {
        //printf("%lu ", global_bf[c]->bytes_per_chip);
        sum += global_bf[c]->bytes_per_chip;
    }
    printf("(SUM: %lu bytes)\n", sum);
    
    sum = 0;
    printf("Sum of bloomfilter assigned bytes per block in one chip: ");
    for(int b=0; b<BLOCK_PER_CHIP; b++) {
        //printf("%lu ", global_bf[0]->bytes_per_blk[b]);
        sum += global_bf[0]->bytes_per_blk[b];
    }
    printf("(SUM: %lu bytes)\n", sum);
    
    sum = 0;
    printf("Sum of bloomfilter assigned bytes per page in one block: ");
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        //printf("%lu ", global_bf[0]->bytes_per_pg[0][p]);
        sum += global_bf[0]->bytes_per_pg[0][p];
    }
    printf("(SUM: %lu bytes)\n", sum);
    
    printf("DONE !!\n");
}

int main() {
    int try_rd = DATA_SET, try_wr = DATA_SET;
    uint32_t lba, pbn;
    double ptv_p=0.0, false_p=0.0;

    read_cnt = write_cnt = true_cnt = false_cnt = found_cnt = notfound_cnt = 0;

    bloom_init();

    /*
     * TRADITIONAL BF
     *
    global_bf = (BF**)malloc(sizeof(BF*) * NUM_BLOCKS);
    for(int i=0; i<NUM_BLOCKS; i++) {
        global_bf[i] = bf_init(NUM_PAGES, 0.1);
    }
    */

    srand((unsigned int)time(NULL));

    printf("TEST DATASET: %d\n", DATA_SET);
    val = data_sz = 0;
    while(try_wr--) {
        lba = generate_lba("WR", "SQ");
        bloom_write(lba, val, "SQ");
        //lba = generate_lba("WR", "RD");
        //bloom_write(lba, val, "RD");
        
        val++; write_cnt++;
    }

    val = 0;
    while(try_rd--) {
        lba = generate_lba("RD", "SQ");
        //lba = generate_lba("RD", "RD");
        bloom_read(lba);

        val++; read_cnt++;
    }

    print_stats();

    bloom_destroy();

    /*
     * TRADITIONAL BF
     *
    for(int i=0; i<NUM_BLOCKS; i++) {
        bf_free(global_bf[i]);
    }
    free(global_bf);
    */

    return 0;
}
