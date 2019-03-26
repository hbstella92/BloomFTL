#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "settings.h"
#include "lsm_settings.h"
#include "bloomfilter.h"
#include "sha256.h"
#include "zlib.h"
#include "zlib_2.h"
#include "zconf.h"

#define CHANNEL 8
#define WAY 8
#define CHIP ((CHANNEL)*(WAY))

#define BLOCK_PER_CHIP 1
#define PAGE_PER_BLOCK 128
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP))
#define TOTAL_PAGE ((TOTAL_BLOCK)*(PAGE_PER_BLOCK))
#define PAGE_SZ (4*(K))
#define CAPACITY ((PAGE_SZ)*(TOTAL_PAGE))
#define DATA_SET (TOTAL_PAGE)
#define PFTL ((PAGE_PER_BLOCK)*(32))

#define LSB (int)(log(TOTAL_BLOCK) / log(2))
#define MASK ((1 << LSB) - 1)

// Compression buffer
#define WRKMEM_SZ 10

// Symbol buffer
#define MAX_SYMBUF 20

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
    BF** bfchip_arr;
} GlobalBF;

typedef struct {
    uint64_t* bits_per_chip;
    uint64_t** bits_per_blk;
    uint64_t*** bits_per_pg;
    uint64_t*** bytes_arr;
} BFManager;

GlobalBF** global_bf;
BFManager* bf_man;

Storage storage;
Block** data_block;

double* fp_arr;

// Compression-related
typedef struct {
    char* comp_arr;
    size_t len;
} Compbuf;

Compbuf*** compbuf;
Compbuf*** decompbuf;
Compbuf*** checkbuf;

typedef struct {
    uint64_t sym_bits_chip; // 1 chip
    uint64_t sym_bits_blk; // 1 block
    uint64_t* sym_bits_pg;
    uint64_t sym_bits_total;
} SManager;

// Symbol (hash func == 1)
char*** symbol;
uint64_t* start;
SManager* st_man;

int* gidx_arr;
int* iidx_arr;
int bitcnt=0;

void bloom_init() {
    double true_p=0.0, false_p=0.0;
   
    // Alloc && Initialize BFManager
    bf_man = (BFManager*)malloc(sizeof(BFManager));
    bf_man->bits_per_chip = (uint64_t*)malloc(sizeof(uint64_t) * CHIP);
    bf_man->bits_per_blk = (uint64_t**)malloc(sizeof(uint64_t*) * CHIP);
    bf_man->bits_per_pg = (uint64_t***)malloc(sizeof(uint64_t**) * CHIP);
    bf_man->bytes_arr = (uint64_t***)malloc(sizeof(uint64_t**) * CHIP);
    for(int c=0; c<CHIP; c++) {
        bf_man->bits_per_blk[c] = (uint64_t*)malloc(sizeof(uint64_t) * BLOCK_PER_CHIP);
        bf_man->bits_per_pg[c] = (uint64_t**)malloc(sizeof(uint64_t*) * BLOCK_PER_CHIP);
        bf_man->bytes_arr[c] = (uint64_t**)malloc(sizeof(uint64_t*) * BLOCK_PER_CHIP);

        bf_man->bits_per_chip[c] = 0;

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            bf_man->bits_per_pg[c][b] = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
            bf_man->bytes_arr[c][b] = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);

            bf_man->bits_per_blk[c][b] = 0;

            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                bf_man->bits_per_pg[c][b][p] = 0;
                bf_man->bytes_arr[c][b][p] = 0;
            }
        }
    }

    // Alloc && Initialize GlobalBF
    global_bf = (GlobalBF**)malloc(sizeof(GlobalBF*) * CHIP);
    for(int c=0; c<CHIP; c++) {
        global_bf[c] = (GlobalBF*)malloc(sizeof(GlobalBF) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_bf[c][b].bfchip_arr = bf_init(1, PAGE_PER_BLOCK);

            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                if(p == 0) {
                    bf_man->bits_per_pg[c][b][p] = 0;
                    bf_man->bytes_arr[c][b][p] = 0;
                }
                else {
                    bf_man->bits_per_pg[c][b][p] = bf_bits(global_bf[c][b].bfchip_arr[p]);

                    int bytes = bf_man->bits_per_pg[c][b][p] / 8;
                    if(bf_man->bits_per_pg[c][b][p] % 8) {
                        bytes++;
                    }
                    bf_man->bytes_arr[c][b][p] = bytes;
                }

                bf_man->bits_per_blk[c][b] += bf_man->bits_per_pg[c][b][p];
            }
            bf_man->bits_per_chip[c] += bf_man->bits_per_blk[c][b];
        }
    }

    // Compression-related
    compbuf = (Compbuf***)malloc(sizeof(Compbuf**) * CHIP);
    decompbuf = (Compbuf***)malloc(sizeof(Compbuf**) * CHIP);
    for(int c=0; c<CHIP; c++) {
        compbuf[c] = (Compbuf**)malloc(sizeof(Compbuf*) * BLOCK_PER_CHIP);
        decompbuf[c] = (Compbuf**)malloc(sizeof(Compbuf*) * BLOCK_PER_CHIP);    
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            compbuf[c][b] = (Compbuf*)malloc(sizeof(Compbuf) * PAGE_PER_BLOCK);
            decompbuf[c][b] = (Compbuf*)malloc(sizeof(Compbuf) * PAGE_PER_BLOCK);
        }
    }

    // Assign different bytes per page-level in all bloomfilters separately
    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                if(p == 0) {
                    compbuf[c][b][p].comp_arr = decompbuf[c][b][p].comp_arr = NULL;
                }
                else {
                    compbuf[c][b][p].comp_arr = (char*)malloc(sizeof(char) * WRKMEM_SZ);
                    decompbuf[c][b][p].comp_arr = (char*)malloc(sizeof(char) * WRKMEM_SZ);
                }
                compbuf[c][b][p].len = decompbuf[c][b][p].len = 0;
            }
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
            free(bf_man->bytes_arr[c][b]);
            free(bf_man->bits_per_pg[c][b]);
            
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                free(compbuf[c][b][p].comp_arr);
                free(decompbuf[c][b][p].comp_arr);
            }

            bf_free(global_bf[c][b].bfchip_arr, PAGE_PER_BLOCK);

            free(compbuf[c][b]);
            free(decompbuf[c][b]);
        }
        free(bf_man->bytes_arr[c]);
        free(bf_man->bits_per_pg[c]);
        free(bf_man->bits_per_blk[c]);
            
        free(global_bf[c]);

        free(compbuf[c]);
        free(decompbuf[c]);
    }
    free(bf_man->bytes_arr);
    free(bf_man->bits_per_pg);
    free(bf_man->bits_per_blk);
    free(bf_man->bits_per_chip);
    free(bf_man);
    free(global_bf);

    free(compbuf);
    free(decompbuf);

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

uint32_t generate_lba(char req, char* pttr) {
    if(req == 'W') {        
        if(!strcmp(pttr, "SEQ")) {
            return val;
        }
        else { // Random and uniform
            uint32_t nodup_rand_val = rand() % DATA_SET;
            int idx;
            while(1) {
                for(idx=0; idx<DATA_SET; idx++) {
                    if(nodup_rand_val == written[idx].lbanum) {
                        nodup_rand_val = rand() % DATA_SET;
                        break;
                    }
                }

                if(nodup_rand_val != written[DATA_SET-1].lbanum) {
                    return nodup_rand_val;
                }
            }
        }
    }
    else if(req == 'R') {
        if(!strcmp(pttr, "SEQ")) {
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

uint64_t cal_combination(uint64_t bits, uint32_t func) {
    uint64_t n = bits; uint64_t r = func;

    if(r > n-r) {
        r = n-r;
    }

    uint64_t up = 1; uint64_t down = 1;
    uint64_t t = r; uint64_t i = 0;
    
    while(t--) {
        up *= (n-i);
        down *= (r-i);
        i++;
    }
    
    return (up / down);
}

void print_byte_as_bit(unsigned char* val, size_t bytes) {
    for(int i=0; i<bytes; i++) {
        for(int j=7; j>=0; j--) {
            printf("%c", (val[i] & (1 << j)) ? '1' : '0');
        } printf(" ");
    } printf("\n");
}

void bloom_write(uint32_t lba, uint32_t value, char* pttr) {
    uint32_t pbn, key, hashkey;
    int chip, way, chnl, blk, offset;
    int global_bf_idx;
    int internal_bf_idx;

    while(1) {
        pbn = (lba & MASK);
        chip = pbn / BLOCK_PER_CHIP;
        way = chip / CHANNEL;
        chnl = chip % CHANNEL;
        blk = pbn % BLOCK_PER_CHIP;
        offset = storage.chip_arr[way][chnl].data_blk[blk].empty;

        if(offset >= PAGE_PER_BLOCK) {
            lba = generate_lba('W', pttr);
        }
        else {
            break;
        }
    }
    
    // Do not have to make bfchip_arr of page 0 (OPTIMIZATION)
    if(offset != 0) {
        key = lba + offset;
        hashkey = hashing_key(key);
        global_bf_idx = bf_set(global_bf[chip][blk].bfchip_arr, offset, hashkey);
        internal_bf_idx = global_bf_idx - global_bf[chip][blk].bfchip_arr[offset]->start;

/*
 * Compression with SymbolTable
 */
        // Encoding BF to Symbol with h
        int symb_end_gidx = start[offset] + st_man->sym_bits_pg[offset] - 1;
        
        int byte, bit;
        while(internal_bf_idx != 0) {
            byte = symb_end_gidx / 8;
            bit = 7 - (symb_end_gidx % 8);
            
            if(internal_bf_idx % 2) {
                symbol[chip][blk][byte] |= (1 << bit);
            }
            symb_end_gidx--;
            internal_bf_idx /= 2;
        }
/**/
    }

/*
 * Compression with zlib
 *
    if(offset != 0) {
        int ret_def;
        ret_def = def(global_bf[chip][blk][offset].bfchip_arr->body, \
                bf_man->bytes_arr[chip][blk][offset], compbuf[chip][blk][offset].comp_arr, \
                &compbuf[chip][blk][offset].len, Z_DEFAULT_COMPRESSION);
        
        if(ret_def != Z_OK) {
            zerr(ret_def);
            exit(1);
        }
    }
*/

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
    int internal_bf_idx = 0;
    
    for(int idx=offset-1; idx>0; idx--) {
        key = lba + idx;
        hashkey = hashing_key(key);
/*
 * Decompression with zlib
 *
        if(idx != 0) {
            int ret_inf;
            ret_inf = inf(compbuf[chip][blk][idx].comp_arr, compbuf[chip][blk][idx].len, \
                    decompbuf[chip][blk][idx].comp_arr, &decompbuf[chip][blk][idx].len);
            
            if(ret_inf != Z_OK) {
                zerr(ret_inf);
                exit(1);
            }
        }
*/

/*
 * Decompression with SymbolTable
 */
        // Decoding Symbol to BF's index
        internal_bf_idx=0;
        
        int symb_front_gidx = start[idx];
        int symb_end_gidx = symb_front_gidx + st_man->sym_bits_pg[idx] - 1;
        
        int byte, bit;
        int j=0;
        while(symb_end_gidx >= symb_front_gidx) {
            byte = symb_end_gidx / 8;
            bit = 7 - (symb_end_gidx % 8);
            
            if(symbol[chip][blk][byte] & (1 << bit)) {
                internal_bf_idx += (int)pow(2, j);
            }

            symb_end_gidx--; j++;
        }
        
        // Bloomfilter checking process (modified with Symbol table)
        //if(bf_check(global_bf[chip][blk].bfchip_arr, idx, hashkey) == true) {
        if(symbol_check(global_bf[chip][blk].bfchip_arr, idx, hashkey, \
                    internal_bf_idx, (int)global_bf[chip][blk].bfchip_arr[idx]->start) == true) {
            if(storage.chip_arr[way][chnl].data_blk[blk].oob[idx] == lba) {
                found_cnt++;
                reading[read_cnt].lbanum = lba;
                reading[read_cnt].level = offset - idx;
                reading[read_cnt].found++;
                return;
            } else {
                notfound_cnt++;
                reading[read_cnt].notfound++;
                continue;
            }
        } else {
            false_cnt++;
        }
/**/

/*
 * Original BF
 *
        // Bloomfilter says true
        if(bf_check(global_bf[chip][blk].bfchip_arr, idx, hashkey) == true) {
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
**/
    }

    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[blk].oob[0] != lba) {
        printf("This should not happen !!\n");
        exit(1);
    } else {
        found_cnt++;
        reading[read_cnt].lbanum = lba;
        reading[read_cnt].level = offset - 0;
        reading[read_cnt].found++;
    }
}

void symbol_init() {
    // Allocate SManager
    st_man = (SManager*)malloc(sizeof(SManager));
    st_man->sym_bits_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    
    memset(st_man->sym_bits_pg, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    st_man->sym_bits_total = st_man->sym_bits_chip = st_man->sym_bits_blk = 0;

    // Allocate symbol delimiter
    start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        start[p] = 0;
    }
    
    // Calculate symbol bits
    uint64_t symbol_bit=0, sum=0;

    st_man->sym_bits_pg[0] = 0;
    start[0] = 0;
    for(int p=1; p<PAGE_PER_BLOCK; p++) {
        symbol_bit = bf_man->bits_per_pg[0][0][p];
        symbol_bit = ceil(log(symbol_bit) / log(2));
        
        st_man->sym_bits_pg[p] = symbol_bit;
        start[p] = sum;
        sum += symbol_bit;
    }
    st_man->sym_bits_blk = sum;
    st_man->sym_bits_chip = sum * BLOCK_PER_CHIP;
    st_man->sym_bits_total = st_man->sym_bits_chip * CHIP;
    
    uint64_t symbol_sum_byte = sum / 8;
    if(sum % 8) {
        symbol_sum_byte++;
    }

    // Allocate symbol table
    symbol = (char***)malloc(sizeof(char**) * CHIP);
    memset(symbol, 0, sizeof(char**) * CHIP);
    for(int c=0; c<CHIP; c++) {
        symbol[c] = (char**)malloc(sizeof(char*) * BLOCK_PER_CHIP);
        memset(symbol[c], 0, sizeof(char*) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            symbol[c][b] = (char*)malloc(sizeof(char) * symbol_sum_byte);
            memset(symbol[c][b], 0, sizeof(char) * symbol_sum_byte);
        }
    }

/*
            // Calculate symbol table size
            int bit_count=0, symbol_bit, sum_of_bits=0;
            uint32_t func, comb;
            BF* prev_bf = global_bf[c][b].bfchip_arr[0];
            uint64_t prev_bits = 0;
            int idx=0;

            for(int p=1; p<PAGE_PER_BLOCK; p++) {
                BF* cur_bf = global_bf[c][b].bfchip_arr[p];
                uint64_t cur_bits = bf_bits(cur_bf);

                if(cur_bits != prev_bits) {
                    if(p == 1) {
                        bit_count = 1;
                        prev_bits = cur_bits;
                        prev_bf = cur_bf;
                        idx = 1;
                        continue;
                    }

                    func = bf_func(prev_bf);
                    comb = cal_combination(prev_bits, func);
                    symbol_bit = ceil(log(comb) / log(2));
                    sum_of_bits += (bit_count * symbol_bit);

                    for(int s=idx, bc=0; s<p; s++, bc++) {
                        sym_table[c][b].start[s] = sum_of_bits - ((bit_count - bc) * symbol_bit);
                    }

                    bit_count = 1;
                    prev_bits = cur_bits;
                    prev_bf = cur_bf;
                    idx = p;
                } else {
                    bit_count++;
                }

                if(p == PAGE_PER_BLOCK-1) {
                    func = bf_func(prev_bf);
                    comb = cal_combination(prev_bits, func);
                    symbol_bit = ceil(log(comb) / log(2));
                    sum_of_bits += (bit_count * symbol_bit);

                    for(int s=idx, bc=0; s<=p; s++, bc++) {
                        sym_table[c][b].start[s] = sum_of_bits - ((bit_count - bc) * symbol_bit);
                    }
                }
            }
        }
    }
*/
}

void symbol_destroy() {
    free(st_man->sym_bits_pg);
    free(st_man);

    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            free(symbol[c][b]);
        }

        free(symbol[c]);
    }
    free(symbol);

    free(start);
}

void print_stats(char* w_type, char* r_type) {
    uint64_t sum = 0;

/*
    printf("\nFOUND LEVEL OF ALL READ REQUESTS\n");
    for(int i=0; i<DATA_SET; i++) {
        printf("req_num: %d level: %d\n", i, reading[i].level);
    }

    printf("\nFALSE POSITIVE RATE OF ALL READ REQUESTS\n");
    for(int i=0; i<DATA_SET; i++) {
        reading[i].fp = (double)reading[i].notfound / (reading[i].found + reading[i].notfound);
        printf("req_num: %d false_p: %.2lf\n", i, reading[i].fp);
    }
*/    
    printf("TEST DATASET: %d\n", DATA_SET);
    printf("TEST TYPE: %s write, %s read\n", w_type, r_type);

    printf("\n### BENCHMARK RESULTS ###\n");
    printf("Total read requests: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("Total false num: %d\n", false_cnt);
    printf("RAF: %.2f\n", (double)(found_cnt + notfound_cnt) / read_cnt);

    uint64_t bytes = 0;
    uint64_t total_bits = 0;
    
    printf("\n### BLOOMFILTER INFO ###\n");
    printf("Sum of BF assigned bits in all chips: ");
    for(int c=0; c<CHIP; c++) {
        sum += bf_man->bits_per_chip[c];
    }
    bytes = sum / 8;
    if(sum % 8) {
        bytes++;
    }
    printf("[%lu bits, %lu bytes]\n", sum, bytes);
    total_bits = sum;
    
    sum = 0;
    printf("Sum of BF assigned bits in all blocks in 1 chip: ");
    for(int b=0; b<BLOCK_PER_CHIP; b++) {
        sum += bf_man->bits_per_blk[0][b];
    }
    bytes = sum / 8;
    if(sum % 8) {
        bytes++;
    }
    printf("[%lu bits, %lu bytes]\n", sum, bytes);
    
    sum = 0;
    printf("Sum of BF assigned bits in all pages in 1 block: ");
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        sum += bf_man->bits_per_pg[0][0][p];
    }
    bytes = sum / 8;
    if(sum % 8) {
        bytes++;
    }
    printf("[%lu bits, %lu bytes] (%.2lf%% of PFTL)\n", sum, bytes, (double)sum/PFTL*100);
    
    bytes = total_bits / 8;
    if(total_bits % 8) {
        bytes++;
    }
    printf("ALL of BF assigned: %lu bits, %lu bytes\n", total_bits, bytes);
    
    printf("\n### SYMBOL TABLE INFO ###\n");
    uint64_t targetsize = st_man->sym_bits_total / 8;
    if(st_man->sym_bits_total % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all chips: [%lu bits, %lu bytes]\n", st_man->sym_bits_total, targetsize);

    targetsize = st_man->sym_bits_chip / 8;
    if(st_man->sym_bits_chip % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all blocks in 1 chip: [%lu bits, %lu bytes]\n", st_man->sym_bits_chip, targetsize);

    targetsize = st_man->sym_bits_blk / 8;
    if(st_man->sym_bits_blk % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all pages in 1 block: [%lu bits, %lu bytes]", st_man->sym_bits_blk, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)st_man->sym_bits_blk/PFTL*100);

    targetsize = st_man->sym_bits_total / 8;
    if(st_man->sym_bits_total % 8) {
        targetsize++;
    }
    printf("ALL of ST assigned: %lu bits, %lu bytes\n", st_man->sym_bits_total, targetsize);
    printf("\nTEST COMPLETE !!\n");
    fflush(stdout);
}

int main(int argc, char** argv) {
    int try_rd = DATA_SET, try_wr = DATA_SET;
    uint32_t lba, pbn;
    double ptv_p=0.0, false_p=0.0;

    read_cnt = write_cnt = true_cnt = false_cnt = found_cnt = notfound_cnt = 0;

    if(argc != 3) {
        printf("Usage: ./simulationFTL W-type R-type\n(Type must be either SEQ or RAND)\n");
        exit(1);
    }

    bloom_init();
    symbol_init();

/*
 * TRADITIONAL BF
 *
    global_bf = (BF**)malloc(sizeof(BF*) * NUM_BLOCKS);
    for(int i=0; i<NUM_BLOCKS; i++) {
        global_bf[i] = bf_init(NUM_PAGES, 0.1);
    }
*/

    printf("TEST START !!\n");
    srand((unsigned int)time(NULL));

    val = data_sz = 0;
    while(try_wr--) {
        lba = generate_lba('W', argv[1]);
        bloom_write(lba, val, argv[1]);
        
        val++; write_cnt++;
    }

    val = 0;
    while(try_rd--) {
        lba = generate_lba('R', argv[2]);
        bloom_read(lba);

        val++; read_cnt++;
    }

    print_stats(argv[1], argv[2]);

    symbol_destroy();
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
