#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "settings.h"
#include "lsm_settings.h"
#include "bloomfilter.h"
#include "sha256.h"
//#include "sha256-arm.h"

#define BLK 0
#define SUPERBLK 1

#define CHANNEL 1
#define WAY 1
#define CHIP ((CHANNEL)*(WAY)) // 1

// block-related
#define BLOCK_PER_CHIP 4
#define PAGE_PER_BLOCK 256 // logical unit
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP)) // 4
#if BLK
#define TOTAL_PAGE ((TOTAL_BLOCK)*(PAGE_PER_BLOCK)) // 1024
#define LP_SZ (4*(K)) // logical page size
#define CAPACITY ((LP_SZ)*(TOTAL_PAGE))
#define PFTL ((PAGE_PER_BLOCK)*32)

#define DATA_SET (TOTAL_PAGE)

/* Buffered write
#define BUF_SZ 4
#define PP_SZ ((LP_SZ)*(BUF_SZ)) // physical page size
#define BLOCK_PPN ((PAGE_PER_BLOCK)/(BUF_SZ)) // 128
*/

#define REBLOOM (int)((PAGE_PER_BLOCK)*0.5)
#endif

// superblock-related
#define SUPERBLK_SZ 4
#define SUPERBLK_PER_CHIP ((BLOCK_PER_CHIP)/(SUPERBLK_SZ)) // 1
#define PAGE_PER_SUPERBLK ((PAGE_PER_BLOCK)*(SUPERBLK_SZ)) // 1024
#define TOTAL_SUPERBLK ((TOTAL_BLOCK)/(SUPERBLK_SZ)) // 1
#if SUPERBLK
#define TOTAL_PAGE ((TOTAL_SUPERBLK)*(PAGE_PER_SUPERBLK)) // 1024
#define LP_SZ (8*(K))
#define CAPACITY ((LP_SZ)*(TOTAL_PAGE))
#define PFTL ((PAGE_PER_SUPERBLK)*32)

#define DATA_SET (TOTAL_PAGE)

#define REBLOOM (int)((PAGE_PER_SUPERBLK)*0.5)
#endif

// Options for lpa generation
#define W_UNIQUE 0x1
#define R_UNIQUE 0x2
#define R_RD_ORDER 0x4
#define R_RD_GEN 0x8

double read_check=0.0;
int read_loop=0;
double w_time=0.0;
double r_time=0.0;

uint32_t val;

uint32_t mask;

int read_cnt;
int write_cnt;
int true_cnt;
int false_cnt;
int found_cnt;
int notfound_cnt;

int disk_write_cnt;

typedef struct {
    int* num_bf; // 0 ~ 256(or 260)
} SBlkManager;

SBlkManager* sblk_man;

/*
typedef struct {
    uint32_t* lpa;
    uint32_t* value;
    int buf_sz;
} BlockBuffer; // TODO: change to struct Page

BlockBuffer** block_buffer;
*/

typedef struct {
    uint32_t page;
    uint32_t oob;
} Page;

typedef struct {
    Page* page_arr;
#if BLK
    int empty; // offset of physical unit
#endif
} Block;

typedef struct {
#if BLK
    Block* data_blk;
#elif SUPERBLK
    Block** data_blk;
    int* empty;
#endif
} Chip;

typedef struct {
    Chip** chip_arr;
} Storage;

typedef struct {
    BF** bfchip_arr;
} GlobalBF;

typedef struct {
    uint64_t* bits_per_pg;
    uint64_t* bytes_arr;
} BFManager;

GlobalBF** global_bf;
BFManager* bf_man;

Storage storage;

// Symbol-related
typedef struct {
    uint8_t*** symbol;

    int** lpa_flag;// TODO: must change to bit array !!!!!
    int** ppa_flag;// TODO: must change to bit array !!!!!
} GlobalSymb;

GlobalSymb global_symb;

typedef struct {
    uint64_t sym_bytes_total;
    uint64_t dead_bytes_total;

    uint64_t sym_bits_chip; // 1 chip
    uint64_t sym_bits_blk; // 1 block
    uint64_t sym_bits_super_blk; // 1 superblock
    uint64_t* sym_bits_pg;
    uint64_t sym_bits_total;
    uint64_t* new_pidx_start;

    uint64_t* sym_start; // bit unit
} SManager;

SManager* st_man;

void bloom_init() {
    double true_p=0.0, false_p=0.0;

#if BLK
    mask = (int)(log(PAGE_PER_BLOCK) / log(2));
#elif SUPERBLK
    mask = (int)(log(PAGE_PER_SUPERBLK) / log(2));
#endif

    // Alloc && Initialize SBlkManager
    sblk_man = (SBlkManager*)malloc(sizeof(SBlkManager));
#if BLK
    sblk_man->num_bf = (int*)malloc(sizeof(int) * TOTAL_BLOCK);
    for(int i=0; i<TOTAL_BLOCK; i++) {
        sblk_man->num_bf[i] = 0;
    }
#elif SUPERBLK
    sblk_man->num_bf = (int*)malloc(sizeof(int) * TOTAL_SUPERBLK);
    for(int i=0; i<TOTAL_SUPERBLK; i++) {
        sblk_man->num_bf[i] = 0;
    }
#endif
   
    // Alloc && Initialize BFManager
    bf_man = (BFManager*)malloc(sizeof(BFManager));
#if BLK
    bf_man->bits_per_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    bf_man->bytes_arr = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
#elif SUPERBLK
    bf_man->bits_per_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SUPERBLK);
    bf_man->bytes_arr = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SUPERBLK);
#endif
   
#if BLK
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
#elif SUPERBLK
    for(int p=0; p<PAGE_PER_SUPERBLK; p++) {
#endif
        bf_man->bits_per_pg[p] = 0;
        bf_man->bytes_arr[p] = 0;
    }

    // Alloc && Initialize global_bf
    global_bf = (GlobalBF**)malloc(sizeof(GlobalBF*) * CHIP);
    for(int c=0; c<CHIP; c++) {
#if BLK
        global_bf[c] = (GlobalBF*)malloc(sizeof(GlobalBF) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_bf[c][b].bfchip_arr = bf_init(1, PAGE_PER_BLOCK);
        }
#elif SUPERBLK
        global_bf[c] = (GlobalBF*)malloc(sizeof(GlobalBF) * SUPERBLK_PER_CHIP);

        for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
            global_bf[c][sb].bfchip_arr = bf_init(1, PAGE_PER_SUPERBLK);
        }
#endif
    }

#if BLK
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
#elif SUPERBLK
    for(int p=0; p<PAGE_PER_SUPERBLK; p++) {
#endif
        if(p != 0) {
            bf_man->bits_per_pg[p] = bf_bits(global_bf[0][0].bfchip_arr[p]);
            
            int targetsize = bf_man->bits_per_pg[p] / 8;
            if(bf_man->bits_per_pg[p] % 8) {
                targetsize++;
            }
            bf_man->bytes_arr[p] = targetsize;
        }
    }

    // Alloc memory for all blocks and pages
    storage.chip_arr = (Chip**)malloc(sizeof(Chip*) * WAY);
    for(int w=0; w<WAY; w++) {
        storage.chip_arr[w] = (Chip*)malloc(sizeof(Chip) * CHANNEL);

        for(int ch=0; ch<CHANNEL; ch++) {
#if BLK
            storage.chip_arr[w][ch].data_blk = (Block*)malloc(sizeof(Block) * BLOCK_PER_CHIP);
#elif SUPERBLK
            storage.chip_arr[w][ch].data_blk = (Block**)malloc(sizeof(Block*) * SUPERBLK_PER_CHIP);
            storage.chip_arr[w][ch].empty = (int*)malloc(sizeof(int) * SUPERBLK_PER_CHIP);
#endif

#if BLK
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                storage.chip_arr[w][ch].data_blk[b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
     
                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    storage.chip_arr[w][ch].data_blk[b].page_arr[p].page = 
                        storage.chip_arr[w][ch].data_blk[b].page_arr[p].oob = 99999;
                }

                storage.chip_arr[w][ch].data_blk[b].empty = 0;
            }
#elif SUPERBLK
            for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
                storage.chip_arr[w][ch].data_blk[sb] = (Block*)malloc(sizeof(Block) * SUPERBLK_SZ);

                storage.chip_arr[w][ch].empty[sb] = 0;

                for(int b=0; b<SUPERBLK_SZ; b++) {
                    storage.chip_arr[w][ch].data_blk[sb][b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

                    for(int p=0; p<PAGE_PER_BLOCK; p++) {
                        storage.chip_arr[w][ch].data_blk[sb][b].page_arr[p].page = 
                        storage.chip_arr[w][ch].data_blk[sb][b].page_arr[p].oob = 99999;
                    }
                }
            }
#endif
        }
    }

    // Alloc memory for block buffer
/*
    block_buffer = (BlockBuffer**)malloc(sizeof(BlockBuffer*) * CHIP);
    for(int c=0; c<CHIP; c++) {
        block_buffer[c] = (BlockBuffer*)malloc(sizeof(BlockBuffer) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            block_buffer[c][b].lpa = (uint32_t*)malloc(sizeof(uint32_t) * BUF_SZ);
            memset(block_buffer[c][b].lpa, 0, sizeof(uint32_t) * BUF_SZ);
            
            block_buffer[c][b].value = (uint32_t*)malloc(sizeof(uint32_t) * BUF_SZ);
            memset(block_buffer[c][b].value, 0, sizeof(uint32_t) * BUF_SZ);

            block_buffer[c][b].buf_sz = 0;
        }
    }
*/
}

void bloom_destroy() {
    free(sblk_man->num_bf);
    free(sblk_man);

    for(int c=0; c<CHIP; c++) {
#if BLK
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            bf_free(global_bf[c][b].bfchip_arr, PAGE_PER_BLOCK);
#elif SUPERBLK
        for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
            bf_free(global_bf[c][sb].bfchip_arr, PAGE_PER_SUPERBLK);
#endif
/*
            free(block_buffer[c][b].lpa);
            free(block_buffer[c][b].value);
*/
        }

        free(global_bf[c]);
/*
        free(block_buffer[c]);
*/
    }

    free(bf_man->bytes_arr);
    free(bf_man->bits_per_pg);
    free(bf_man);
    free(global_bf);
/*
    free(block_buffer);
*/

    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
#if BLK
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                free(storage.chip_arr[w][ch].data_blk[b].page_arr);
#elif SUPERBLK
            for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
                for(int b=0; b<SUPERBLK_SZ; b++) {
                    free(storage.chip_arr[w][ch].data_blk[sb][b].page_arr);
                }

                free(storage.chip_arr[w][ch].data_blk[sb]);
#endif
            }
            
            free(storage.chip_arr[w][ch].data_blk);
        }

        free(storage.chip_arr[w]);
    }

    free(storage.chip_arr);
}

/* SHA256 by hanbyeol
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
*/

/* SHA256 by jiho
uint32_t hashing_key(uint32_t key) {
	uint32_t hashkey; 
	uint32_t state[8];
	
	sha256_calculate(state, key);

	hashkey = state[0];
	for(int i=1; i<8; i++) {
		hashkey ^= state[i];
	}

	return hashkey;
}
*/

// Fibonacci hash
static inline uint32_t hashing_key(uint32_t key) {
    return (uint32_t)((0.618033887 * key) * 1024);
}

int compare(const void *a, const void *b) {
	uint32_t num1 = *(uint32_t *)a;
	uint32_t num2 = *(uint32_t *)b;

	if (num1 < num2)
		return -1;
	if (num1 > num2)
		return 1;
	return 0;
}

void shuffle(uint32_t *arr, int length) {
	int i, j, t;
	for(i=0; i<length-1; i++) {
		j = i + rand() / (RAND_MAX / (length - i) + 1);
		t = arr[j];
		arr[j] = arr[i];
		arr[i] = t;
	}
}

/* 
 *
    NOTE: If SEQ_FLAG is delivered, these flags does not work!
        W_UNIQUE: Random && Unique
        R_UNIQUE: Random && Unique (If W_UNIQUE is not set, R_UNIQUE doesn't work)
        R_RD_ORDER:
            ex) If W_UNIQUE is set and writes lpa 3-1-2, R_RD_ORDER reads lpa 2-1-3
        R_RD_GEN: Random && Duplicated
 */
void make_test_set(uint32_t *w_arr, uint32_t *r_arr, char *w_t, char *r_t, uint8_t options, uint32_t test_size) {
    uint8_t *page_usage;
    uint32_t make_cnt, lpa, pbn;
    uint8_t read = 1, write = 1;

    if(!strcmp(w_t, "RAND")) {
        write = 0;
    }

    // Generate write lpa
    if(!write && !(options & W_UNIQUE)) {
        // Check whether the physical block corresponding to lpa is full or not (Because there is no GC algorithm)
        make_cnt = 0;
        page_usage = (uint8_t*)malloc(TOTAL_PAGE);
        memset(page_usage, 0, TOTAL_PAGE);

        while(make_cnt < test_size) {
            lpa = rand() % TOTAL_PAGE;
            pbn = lpa >> mask; // random write

            if(page_usage[pbn] == PAGE_PER_BLOCK) {
                continue;
            }

            w_arr[make_cnt++] = lpa;
            page_usage[pbn]++;
        }

        free(page_usage);
        /*
        // After GC is implemented
        for(int i = 0; i < test_size; i++) {
        w_arr[i] = rand() % TOTAL_PAGE;
        }
         */
    }
    else {
        for(int i = 0; i < test_size; i++) {
            w_arr[i] = i;
        }

        if(!write) {
            shuffle(w_arr, test_size);
        }
    }
    if(!strcmp(r_t, "RAND")) {
        read = 0;
    }

    // Generate read lpa
    if(read) { // Seq read
        memcpy(r_arr, w_arr, sizeof(uint32_t) * test_size);
        if(!write) {
            qsort(r_arr, test_size, sizeof(uint32_t), compare);
        }
    }
    else { // Rand read
        if((write && !(options & R_UNIQUE)) ||
                (!write && ((options & (W_UNIQUE | R_UNIQUE)) == W_UNIQUE ||\
                            (options & (W_UNIQUE | R_RD_GEN)) == R_RD_GEN))) {

            for(int i = 0; i < test_size; i++) {
                r_arr[i] = w_arr[rand() % test_size];
            }
        }
        else {
            memcpy(r_arr, w_arr, sizeof(uint32_t) * test_size);

            if(write || (options & R_RD_ORDER)) {
                shuffle(r_arr, test_size);
            }
        }
    }
}

// For debugging
void print_byte_as_bit(unsigned char* val, size_t bytes) {
    for(int i=0; i<bytes; i++) {
        for(int j=7; j>=0; j--) {
            printf("%c", (val[i] & (1 << j)) ? '1' : '0');
        } printf(" ");
    } printf("\n");
}

int lpa_compare(const void* a, const void* b) {
    Page* one = (Page*)a;
    Page* two = (Page*)b;

    if(one->oob < two->oob) {
        return -1;
    }
    else if(one->oob > two->oob) {
        return 1;
    }
    else {
        return 0;
    }
}

#if BLK
void symbol_resymbolize(uint32_t pbn, int chip, int way, int chnl, int blk, int valid_start, int num_flush) {
#elif SUPERBLK
void symbol_resymbolize(uint32_t pbn, int chip, int way, int chnl, int blk, int valid_start, int num_flush, int sb) {
#endif
    uint32_t hashkey;
#if BLK
    storage.chip_arr[way][chnl].data_blk[blk].empty += num_flush;
    
    free(global_symb.symbol[chip][blk]);
    sblk_man->num_bf[pbn] = 0;

    global_symb.symbol[chip][blk] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
    memset(global_symb.symbol[chip][blk], 0, sizeof(uint8_t) * st_man->sym_bytes_total);

    int ppa_end_idx = storage.chip_arr[way][chnl].data_blk[blk].empty;
    int new_bf_idx = 1;
    
    // Set new symbol and Set dead bit for lpa and ppa (TODO: BF's new start idx is 1)
    for(int i=0; i<ppa_end_idx; i++) {
        if(i < valid_start) {
            if(global_symb.ppa_flag[pbn][i] == 1) {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob, new_bf_idx);
                new_bf_idx++;
            }
        }
        else if(i > valid_start) {
            if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i-1].oob + 1 ==
                    storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob) {
                global_symb.ppa_flag[pbn][i] = 0;
            }
            else {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob, new_bf_idx);

                global_symb.ppa_flag[pbn][i] = 1;

                new_bf_idx++;
            }
        }
        else {
            if(i == 0) {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob, new_bf_idx);

                global_symb.ppa_flag[pbn][i] = 1;

                new_bf_idx++;
            }
            else {
                if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i-1].oob + 1 ==
                        storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob) {
                    global_symb.ppa_flag[pbn][i] = 0;
                }
                else {
                    hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                    symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                            global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob, new_bf_idx);

                    global_symb.ppa_flag[pbn][i] = 1;

                    new_bf_idx++;
                }
            }
        }
    }

    sblk_man->num_bf[pbn] = new_bf_idx - 1;
#elif SUPERBLK
    int superblk = chip * SUPERBLK_PER_CHIP + sb;

    storage.chip_arr[way][chnl].empty[sb] += num_flush;

    free(global_symb.symbol[chip][sb]);
    sblk_man->num_bf[superblk] = 0;

    global_symb.symbol[chip][sb] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
    memset(global_symb.symbol[chip][sb], 0, sizeof(uint8_t) * st_man->sym_bytes_total);

    int ppa_end_idx = storage.chip_arr[way][chnl].empty[sb];
    int new_bf_idx = 1;

    // Set new symbol and Set dead bit for lpa and ppa (TODO: BF's new start idx is 1)
    for(int pa=0, b=0; pa<ppa_end_idx, b<SUPERBLK_SZ; pa++) {
        bool make_new_bf = false;

        int local_pa = pa % PAGE_PER_BLOCK;

        if(pa < valid_start) {
            if(global_symb.ppa_flag[superblk][pa] == 1) {
                make_new_bf = true;
            }
        }
        else if(pa > valid_start) {
            if(local_pa == 0) {
                if(storage.chip_arr[way][chnl].data_blk[sb][b-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
                        storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[0].oob) {
                    global_symb.ppa_flag[superblk][pa] = 0;
                }
                else {
                    make_new_bf = true;
                }
            }
            else {
                if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 ==
                        storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
                    global_symb.ppa_flag[superblk][pa] = 0;
                }
                else {
                    make_new_bf = true;
                }
            }
        }
        else {
            if(pa == 0) {
                make_new_bf = true;
            }
            else {
                if(global_symb.ppa_flag[superblk][pa-1] == 1) {
                    if(local_pa == 0) {
                        if(storage.chip_arr[way][chnl].data_blk[sb][b-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
                                storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[0].oob) {
                            global_symb.ppa_flag[superblk][pa] = 0;
                        }
                        else {
                            make_new_bf = true;
                        }
                    }
                    else {
                        if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 == 
                                storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
                            global_symb.ppa_flag[superblk][pa] = 0;
                        }
                        else {
                            make_new_bf = true;
                        }
                    }
                }
                else {
                    make_new_bf = true;
                }
            }
        }

        if(make_new_bf == true) {
            hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob);
            symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                    global_symb.symbol[chip][sb], st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob, new_bf_idx);

            global_symb.ppa_flag[superblk][pa] = 1;

            new_bf_idx++;
        }

        if(!((pa + 1) % PAGE_PER_BLOCK)) {
            b++;
        }

        if(pa + 1 == storage.chip_arr[way][chnl].empty[sb]) {
            break;
        }
    }
    
    sblk_man->num_bf[superblk] = new_bf_idx - 1;
#endif
}

#if BLK
void bloom_gc(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, Page* evict_list, int num_evict) {
#elif SUPERBLK
void bloom_gc(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, Page* evict_list, int num_evict, int* sb) {
    Page *candi_list;
    int num_candi = 0;
#endif
    Page *valid_list, prev_page;
    uint32_t hashkey;
    int num_valid = 0;
#if BLK
    valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

    // Copy valid pages(oob related) in evict block
    uint32_t lpa_start = (*pbn) * PAGE_PER_BLOCK;
    for(int i=0; i<PAGE_PER_BLOCK; i++) {
        if(global_symb.lpa_flag[*pbn][i] == 1) {
            valid_list[num_valid++].oob = lpa_start + i;
        }
    }

    int end_idx = storage.chip_arr[*way][*chnl].data_blk[*blk].empty;
    for(int i=0; i<num_valid; i++) {
        for(int j=end_idx-1; j>=0; j--) {
            if(valid_list[i].oob == storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[j].oob) {
                valid_list[i].page = storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[j].page;
                break;
            }
        }
    }

    // Free exist block
    free(storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr);

    // Re-allocate and Initialize new block
    storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p].page = 
            storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p].oob = 99999;
    }
    storage.chip_arr[*way][*chnl].data_blk[*blk].empty = 0;
    
    memset(global_symb.ppa_flag[*pbn], 0, sizeof(int) * PAGE_PER_BLOCK);
    
    // Flush valid pages
    memcpy(storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr, valid_list, sizeof(Page) * num_valid);
    symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, 0, num_valid);
    
    disk_write_cnt += num_valid;
#elif SUPERBLK
    if(!evict_list) {
        int superblk = (*chip) * SUPERBLK_PER_CHIP + (*sb);

        candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
        valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

        memset(candi_list, 0, sizeof(Page) * PAGE_PER_BLOCK);
        memset(valid_list, 0, sizeof(Page) * PAGE_PER_BLOCK);

        // Copy valid pages in evict block
        for(int i=PAGE_PER_BLOCK-1; i>=0; i--) {
            bool is_valid = true;

            for(int b=SUPERBLK_SZ-1; b>0; b--) {
                for(int p=PAGE_PER_BLOCK-1; p>=0; p--) {
                    if(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i].oob == 
                            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                        is_valid = false;
                        break;
                    }
                }

                if(!is_valid) {
                    break;
                }
            }

            if(is_valid == true) {
                candi_list[num_candi++] = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i];
            }
        }

        qsort(candi_list, num_candi, sizeof(Page), lpa_compare);

        // Remove duplicate
        prev_page = valid_list[0] = candi_list[0];
        num_valid++;
        for(int i=1; i<num_candi; i++) {
            if(prev_page.oob != candi_list[i].oob) {
                valid_list[num_valid++] = candi_list[i];
            }

            prev_page = candi_list[i];
        }

        // Free exist block
        free(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr);

        // Re-allocate and Initialize new block
        Page* tmp_page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

        storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[p].page = 
                storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[p].oob = 99999;
        }
        storage.chip_arr[*way][*chnl].empty[*sb] = (3 * PAGE_PER_BLOCK);
        
        int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SUPERBLK);
        memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SUPERBLK);

        memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK], sizeof(int) * PAGE_PER_BLOCK * 3);
        memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SUPERBLK);

        tmp_page_arr = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;
        for(int b=0; b<SUPERBLK_SZ-1; b++) {
            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
                storage.chip_arr[*way][*chnl].data_blk[*sb][b+1].page_arr;
        }
        storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr = tmp_page_arr;

        // Flush valid pages
        memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr, valid_list, sizeof(Page) * num_valid);
        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, PAGE_PER_BLOCK * 3, num_valid, *sb);

        disk_write_cnt += num_valid;

        //free(tmp_page_arr);
        //free(tmp_ppa_flag);
        //free(candi_list);
    }
    else { // During RB, GC triggererd
        int superblk = (*chip) * SUPERBLK_PER_CHIP + (*sb);
        int evict_blk = -1;
        int new_blk_start;

        candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SUPERBLK * 2);
        valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SUPERBLK);

        memset(candi_list, 0, sizeof(Page) * PAGE_PER_SUPERBLK * 2);
        memset(valid_list, 0, sizeof(Page) * PAGE_PER_SUPERBLK);

        // Copy valid pages in evict block
        memcpy(candi_list, evict_list, sizeof(Page) * num_evict);
        //num_candi += num_evict;
// TODO
int candi_sz = num_candi;
        do {
            num_valid = 0;
            memset(valid_list, 0, sizeof(Page) * PAGE_PER_SUPERBLK);

            evict_blk++;

            for(int i=PAGE_PER_BLOCK-1; i>=0; i--) {
                bool is_valid = true;

                for(int b=SUPERBLK_SZ-1; b>evict_blk; b--) {
                    for(int p=PAGE_PER_BLOCK-1; p>=0; p--) {
                        if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i].oob ==
                                storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                            is_valid = false;
                            break;
                        }
                    }

                    if(!is_valid) {
                        break;
                    }
                }

                if(is_valid == true) {
                    candi_list[num_candi++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i];
                }
            }

            qsort(candi_list, num_candi, sizeof(Page), lpa_compare);

            prev_page = valid_list[0] = candi_list[0];
            num_valid++;
            for(int i=1; i<num_candi; i++) {
                if(prev_page.oob != candi_list[i].oob) {
                    valid_list[num_valid++] = candi_list[i];
                }

                prev_page = candi_list[i];
            }

            candi_list = valid_list;
        } while((((SUPERBLK_SZ - 1 - evict_blk) * PAGE_PER_BLOCK) + num_candi > PAGE_PER_SUPERBLK) && (evict_blk + 1 < SUPERBLK_SZ));
/*
        do {
            evict_blk++;

            for(int i=PAGE_PER_BLOCK-1; i>=0; i--) {
                bool is_valid = true;

                for(int b=SUPERBLK_SZ-1; b>evict_blk; b--) {
                    for(int p=PAGE_PER_BLOCK-1; p>=0; p--) {
                        if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i].oob ==
                                storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                            is_valid = false;
                            break;
                        }
                    }

                    if(!is_valid) {
                        break;
                    }
                }

                if(is_valid == true) {
                    candi_list[num_candi++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i];
                }
            }

            //
        } while((((SUPERBLK_SZ - 1 - evict_blk) * PAGE_PER_BLOCK) + num_candi > PAGE_PER_SUPERBLK) && (evict_blk + 1 < SUPERBLK_SZ));

        qsort(candi_list, num_candi, sizeof(Page), lpa_compare);

printf("num_candi: %d\n", num_candi);
for(int i=0; i<num_candi; i++)
    printf("idx %d\tcandi_list: %u\n", i, candi_list[i].oob);
fflush(stdout);

        // Remove duplicate
        prev_page = valid_list[0] = candi_list[0];
        num_valid++;
        for(int i=1; i<num_candi; i++) {
            if(prev_page.oob != candi_list[i].oob) {
                valid_list[num_valid++] = candi_list[i];
            }

            prev_page = candi_list[i];
        }
*/
printf("num_valid: %d\n", num_valid);
for(int i=0; i<num_valid; i++)
    printf("idx: %d\tvalid_list: %u\n", i, valid_list[i].oob);
printf("\n");
fflush(stdout);

        // Free exist block
        evict_blk++;
        new_blk_start = SUPERBLK_SZ - evict_blk;

        for(int b=0; b<evict_blk; b++) {
            free(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr);

            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[p].page = 
                    storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[p].oob = 99999;
            }
        }

printf("before empty: %d\n", storage.chip_arr[*way][*chnl].empty[*sb]);

        storage.chip_arr[*way][*chnl].empty[*sb] -= (PAGE_PER_BLOCK * evict_blk);

printf("after empty: %d\n", storage.chip_arr[*way][*chnl].empty[*sb]);

        // Re-allocate and Initialize new block
        Page* tmp_page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
   
        int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SUPERBLK);
        memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SUPERBLK);

        memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK * evict_blk], \
                sizeof(int) * PAGE_PER_BLOCK * new_blk_start);
        memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SUPERBLK);

        for(int i=0; i<evict_blk; i++) {
            tmp_page_arr = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;
            for(int b=0; b<SUPERBLK_SZ-1; b++) {
                storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
                    storage.chip_arr[*way][*chnl].data_blk[*sb][b + 1].page_arr;
            }
            storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr = tmp_page_arr;
        }

printf("new_blk_start: %d\tevict_blk: %d\n", new_blk_start, evict_blk); fflush(stdout);

        // Flush valid pages
        int remain_valid = num_valid;
        int j = 0; int v = 0;
        int new_ppn = storage.chip_arr[*way][*chnl].empty[*sb] % PAGE_PER_BLOCK;
        int append_blk = storage.chip_arr[*way][*chnl].empty[*sb] / PAGE_PER_BLOCK;

printf("new_ppn: %d\tappend_blk: %d\n\n", new_ppn, append_blk);
fflush(stdout);

        while(remain_valid != 0) {
            if(j == 0) {

printf("append_blk: %d\tv: %d (remain %d)\n", append_blk, 0, remain_valid);
fflush(stdout);
                memcpy(&storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk].page_arr[new_ppn], valid_list, sizeof(Page) * (PAGE_PER_BLOCK - new_ppn));

                remain_valid -= (PAGE_PER_BLOCK - new_ppn);
                v += (PAGE_PER_BLOCK - new_ppn);
            }
            else {
                if(remain_valid > PAGE_PER_BLOCK) {

printf("append_blk + j: %d\tv: %d (remain %d)\n", append_blk + j, v, remain_valid);
fflush(stdout);
                    memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk + j].page_arr, &valid_list[v], sizeof(Page) * PAGE_PER_BLOCK);

                    remain_valid -= PAGE_PER_BLOCK;
                    v += PAGE_PER_BLOCK;
                }
                else {

printf("append_blk + j: %d\tv: %d (remain %d)\n", append_blk + j, v, remain_valid);
fflush(stdout);
                    memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk + j].page_arr, &valid_list[v], sizeof(Page) * remain_valid);

                    remain_valid = 0;
                }
            }

            j++;
        }
        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, PAGE_PER_BLOCK * append_blk + new_ppn, num_valid, *sb);

        disk_write_cnt += num_valid;

        //free(tmp_page_arr);
        //free(tmp_ppa_flag);
        //free(candi_list);
    }

#endif

    free(valid_list);
}

void debugging(int way, int chnl, int blk, int pbn);

#if BLK
void bloom_rebloom(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn) {
#elif SUPERBLK
void bloom_rebloom(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, int* sb) {
#endif
    Page *candi_list, *evict_list, prev_page;
    uint32_t hashkey;
    int num_candi = 0, num_evict = 0, candi_sz = 0;
#if BLK
    num_candi = storage.chip_arr[*way][*chnl].data_blk[*blk].empty;
    
    candi_list = (Page*)malloc(sizeof(Page) * num_candi);
    evict_list = (Page*)malloc(sizeof(Page) * num_candi);

    // Copy valid pages and Sort with lpa order
    for(int i=0; i<num_candi; i++) {
        if(global_symb.ppa_flag[*pbn][i] == 1) {
            candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[i];
        }
    }
    
    qsort(candi_list, candi_sz, sizeof(Page), lpa_compare);

    // Elect rebloomed-pages
    prev_page = candi_list[0];
    for(int i=1; i<candi_sz; i++) {
        if(prev_page.oob == candi_list[i].oob) {
            if(num_evict == 0) {
                evict_list[num_evict++] = candi_list[i];
            }
            else {
                if(evict_list[num_evict - 1].oob == candi_list[i].oob) {
                    evict_list[num_evict - 1] = candi_list[i];
                }
                else {
                    evict_list[num_evict++] = candi_list[i];
                }
            }
        }
        else if(prev_page.oob + 1 == candi_list[i].oob) {
            if(num_evict == 0) {
                evict_list[num_evict++] = prev_page;
                evict_list[num_evict++] = candi_list[i];
            }
            else {
                if(evict_list[num_evict - 1].oob == prev_page.oob) {
                    evict_list[num_evict - 1] = prev_page;
                    evict_list[num_evict++] = candi_list[i];
                }
                else {
                    evict_list[num_evict++] = prev_page;
                    evict_list[num_evict++] = candi_list[i];
                }
            }
        }

        prev_page = candi_list[i];
    }

    // Set invalid bit with elected rebloomed-pages
    for(int i=0; i<num_evict; i++) {
        for(int j=0; j<num_candi; j++) {
            if(storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[j].oob == evict_list[i].oob) {
                global_symb.ppa_flag[*pbn][j] = 0;
            }
        }
    }

    if(num_candi + num_evict >= PAGE_PER_BLOCK) {
        bloom_gc(pbn, chip, way, chnl, blk, ppn, evict_list, num_evict);
        free(evict_list);
        free(candi_list);
        return;
    }

    // Append rebloomed pages
    int remain_evict = num_evict;
    for(int i=0, p=num_candi; i<num_evict; i++, p++) {
        if(p == PAGE_PER_BLOCK) {
            symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, num_evict);

            bloom_gc(pbn, chip, way, chnl, blk, ppn, NULL, 0);

            p = 0;
            num_candi = 0;
            remain_evict = num_evict - i;
        }

        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p] = evict_list[i];

        disk_write_cnt++;
    }

    symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, remain_evict);
#elif SUPERBLK
    int superblk = (*chip) * SUPERBLK_PER_CHIP + (*sb);
    num_candi += storage.chip_arr[*way][*chnl].empty[*sb];

    candi_list = (Page*)malloc(sizeof(Page) * num_candi);
    evict_list = (Page*)malloc(sizeof(Page) * num_candi);

printf("num_candi: %d\n", num_candi); fflush(stdout);

    // Copy valid pages and Sort with lpa order
    for(int b=0; b<=(*blk); b++) {
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            int gidx = b * PAGE_PER_BLOCK + p;
            
            if(global_symb.ppa_flag[superblk][gidx] == 1) {
                candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
            }

            if(gidx == num_candi - 1) {
            //if((gidx == num_candi - 1) && (b == (*blk))) {
                break;
            }
        }
    }

    qsort(candi_list, candi_sz, sizeof(Page), lpa_compare);

    // Elect rebloomed-pages
    prev_page = candi_list[0];
    for(int i=1; i<candi_sz; i++) {
        if(prev_page.oob == candi_list[i].oob) {
            if(num_evict == 0) {
                evict_list[num_evict++] = candi_list[i];
            }
            else {
                if(evict_list[num_evict - 1].oob == candi_list[i].oob) {
                    evict_list[num_evict - 1] = candi_list[i];
                }
                else {
                    evict_list[num_evict++] = candi_list[i];
                }
            }
        }
        else if(prev_page.oob + 1 == candi_list[i].oob) {
            if(num_evict == 0) {
                evict_list[num_evict++] = prev_page;
                evict_list[num_evict++] = candi_list[i];
            }
            else {
                if(evict_list[num_evict - 1].oob == prev_page.oob) {
                    evict_list[num_evict - 1] = prev_page;
                    evict_list[num_evict++] = candi_list[i];
                }
                else {
                    evict_list[num_evict++] = prev_page;
                    evict_list[num_evict++] = candi_list[i];
                }
            }
        }

        prev_page = candi_list[i];
    }

    // Set invalid bit with elected rebloomed-pages
    for(int i=0; i<num_evict; i++) {
        for(int b=0; b<(*blk); b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                int gidx = b * PAGE_PER_BLOCK + p;
                
                if(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob == evict_list[i].oob) {
                    global_symb.ppa_flag[superblk][gidx] = 0;
                }

                // TODO : fix it (think about)
                if(gidx == num_candi - 1) {
                    break;
                }
            }
        }
    }

    if(num_candi + num_evict >= PAGE_PER_SUPERBLK) {
printf("During RB, GC triggered !!\n");
printf("\nnum_evict: %d\n", num_evict);
        
        bloom_gc(pbn, chip, way, chnl, blk, ppn, evict_list, num_evict, sb);

        free(evict_list);
        free(candi_list);
        return;
    }
    
    // Append rebloomed pages
    int remain_evict = num_evict;
    for(int i=0, p=num_candi%PAGE_PER_BLOCK, b=(*blk); i<num_evict; i++, p++) {
        if(p == PAGE_PER_BLOCK) {
            b++;
            p = 0;
        }
        
        storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p] = evict_list[i];
        disk_write_cnt++;
    }

    // TODO: sync with prototype
    symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, remain_evict, *sb);
#endif
    free(evict_list);
    free(candi_list);
}

void debugging(int way, int chnl, int blk, int pbn) {
    int sum_num_bf = 0;
    uint64_t bfsum = 0, targetsize;

    printf("ppn\tppa\twrite_lpa\tppa_flag(BF exists)\tlpa_flag(lpa exists)\tlpa_list\n");
#if BLK
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        printf("%d\t%d\t\t%u\t\t\t", p, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
        printf("%d\t\t\t\t\t%d\t\t\t%u\n", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
        if(global_symb.ppa_flag[pbn][p] == 1) {
            sum_num_bf++;
        }
    }
    printf("result: # of BF: %d\n", sum_num_bf);
    printf("\n");
#elif SUPERBLK
    int chip = way * CHANNEL + chnl;

    for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
        int superblk = chip * SUPERBLK_PER_CHIP + sb;
        sum_num_bf = 0;
      
        printf("[SUPERBLOCK %d in chip %d]\n", sb, way*CHANNEL+chnl);

        for(int b=0; b<SUPERBLK_SZ; b++) {
            printf("[Block %d]\n", b);
            printf("idx\twrite_lpa\tppa_flag(BF)\tlpa_flag\tlpalist\t\tnew_bf_idx\n");

            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                int pa = b * PAGE_PER_BLOCK + p;
                
                printf("%d\t\t%u\t\t\t", p, storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[p].oob);
                printf("%d\t\t\t\t%u\t\t%u\t\t\t", global_symb.ppa_flag[superblk][pa], \
                        global_symb.lpa_flag[superblk][pa], pa);

                if(global_symb.ppa_flag[superblk][pa] == 1) {
                    sum_num_bf++;
                    printf("%d\n", sum_num_bf);
                }
                else {
                    printf("\n");
                }
            } printf("\n");
        } printf("\n");

        printf("Superblock %d (# of valid BFs: %d)\n", sb, sblk_man->num_bf[superblk]);
    } printf("\n");
#endif

#if BLK
    for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
        bfsum += st_man->sym_bits_pg[p];
    }
    bfsum += (2 * PAGE_PER_BLOCK);
    targetsize = bfsum / 8;
    if(bfsum % 8) {
        targetsize++;
    }
    printf("Sum of ST bits in block %u: [%lu bits, %lu bytes]", pbn, bfsum, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)bfsum/PFTL*100);
    printf("\n");
#elif SUPERBLK
    //
#endif
}

// 8K no-buffered write with Reblooming
void bloom_write(uint32_t lpa, uint32_t value, char* pttr) {
    uint32_t pbn, hashkey;
#if SUPERBLK
    uint32_t superblk;
#endif
    int chip, way, chnl, blk, ppn;
    int lpa_start, new_bf_idx;
    Page prev_page;

#if BLK
    pbn = lpa >> mask; // random write
    blk = pbn % BLOCK_PER_CHIP;
    
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    
    ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

    // GC triggered
    if(ppn >= PAGE_PER_BLOCK) {
        printf("\nBefore GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        debugging(way, chnl, blk, pbn);

        bloom_gc(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0);
        ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

        printf("\nAfter GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        debugging(way, chnl, blk, pbn);
    }
    
    // RB triggered (# of valid BFs in one superblock exceeds 50%)
    if(sblk_man->num_bf[pbn] >= REBLOOM) {
        printf("\nBefore RB in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        debugging(way, chnl, blk, pbn);

        bloom_rebloom(&pbn, &chip, &way, &chnl, &blk, &ppn);
        ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

        printf("\nAfter RB in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        debugging(way, chnl, blk, pbn);
    }

    lpa_start = ppn;
    new_bf_idx = sblk_man->num_bf[pbn] + 1;
    bool create_new_bf = false;

    storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start].oob = lpa;
    storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start].page = value;
    write_cnt++;
    disk_write_cnt++;

    global_symb.lpa_flag[pbn][lpa % PAGE_PER_BLOCK] = 1;

    if(lpa_start != 0) {
        // TODO: how to synchronize BF's start? 0 or 1?
        if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start - 1].oob + 1 == lpa) {
            global_symb.ppa_flag[pbn][lpa_start] = 0;
        }
        else {
            create_new_bf = true;
        }
    }

    if(create_new_bf == true) {
        hashkey = hashing_key(lpa);
        symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, \
                hashkey + new_bf_idx, global_symb.symbol[chip][blk], \
                st_man->sym_start, st_man->sym_bits_pg, false);

        sblk_man->num_bf[pbn]++;
        global_symb.ppa_flag[pbn][lpa_start] = 1;
    }

    storage.chip_arr[way][chnl].data_blk[blk].empty++;
#elif SUPERBLK
    int sb;

    superblk = lpa >> mask; // global superblk
    sb = superblk % 2;

    chip = superblk / SUPERBLK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;

    ppn = storage.chip_arr[way][chnl].empty[sb]; // TODO: global ppn
    blk = ppn / PAGE_PER_BLOCK; // in one superblock
    pbn = sb * SUPERBLK_SZ + blk; // in one chip

    // GC triggered
    if(ppn >= PAGE_PER_SUPERBLK) {
        printf("\nBefore GC in block %u of superblock %u (# of valid BFs: %d)\n", blk-1, \
                superblk, sblk_man->num_bf[superblk]);
        debugging(way, chnl, blk, pbn);

printf("Before GC ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);
        bloom_gc(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0, &sb);
        
        ppn = storage.chip_arr[way][chnl].empty[sb];
        blk = ppn / PAGE_PER_BLOCK;
        pbn = sb * SUPERBLK_SZ + blk;
printf("After GC ppn: %d\tblk: %d\tpbn: %d\n\n", ppn, blk, pbn); fflush(stdout);

        printf("\nAfter GC in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
                superblk, sblk_man->num_bf[superblk]);
        debugging(way, chnl, blk, pbn);
    }

    // RB triggered (# of valid BFs in one superblock exceeds 50%)
    if(sblk_man->num_bf[superblk] >= REBLOOM) {
        printf("\nBefore RB in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
                superblk, sblk_man->num_bf[superblk]);
        debugging(way, chnl, blk, pbn);

printf("Before RB ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);
        bloom_rebloom(&pbn, &chip, &way, &chnl, &blk, &ppn, &sb);
        
        ppn = storage.chip_arr[way][chnl].empty[sb];
        blk = ppn / PAGE_PER_BLOCK;
        pbn = sb * SUPERBLK_SZ + blk;
printf("After RB ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);
        
        printf("\nAfter RB in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
                superblk, sblk_man->num_bf[superblk]);
        debugging(way, chnl, blk, pbn);
    }

    int global_lpa_start = ppn;
    lpa_start = global_lpa_start % PAGE_PER_BLOCK;
    new_bf_idx = sblk_man->num_bf[superblk] + 1;
    bool create_new_bf = false;

    storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lpa_start].oob = lpa;
    storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lpa_start].page = value;
    write_cnt++;
    disk_write_cnt++;
    
    global_symb.lpa_flag[superblk][lpa % PAGE_PER_SUPERBLK] = 1;

    if(global_lpa_start != 0) {
        // TODO: how to synchronize BF's start? 0 or 1?
        if(lpa_start == 0) {
            if(storage.chip_arr[way][chnl].data_blk[sb][blk-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 == lpa) {
                global_symb.ppa_flag[superblk][global_lpa_start] = 0;

            }
            else {
                create_new_bf = true;
            }
        }
        else {
            if(storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lpa_start - 1].oob + 1 == lpa) {

                global_symb.ppa_flag[superblk][global_lpa_start] = 0;
            }
            else {
                create_new_bf = true;
            }
        }
    }
    else {
        create_new_bf = true;
    }

    if(create_new_bf == true) {
        hashkey = hashing_key(lpa);
        symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, \
                hashkey + new_bf_idx, global_symb.symbol[chip][sb], \
                st_man->sym_start, st_man->sym_bits_pg, false);

        sblk_man->num_bf[superblk]++;
        global_symb.ppa_flag[superblk][global_lpa_start] = 1;
    }

    storage.chip_arr[way][chnl].empty[sb]++;
#endif
}

// 8K no-buffered write with Reblooming
void bloom_read(uint32_t lpa) {
    uint32_t pbn, hashkey;
#if SUPERBLK
    uint32_t superblk;
#endif
    int chip, way, chnl, blk, ppn;
    struct timeval strt, end;
    int max_seq_cnt = 0;
   
#if BLK
    pbn = lpa >> mask;
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;

    ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

    int new_bf_idx = sblk_man->num_bf[pbn];
    bool rb_or_gc = (new_bf_idx == ppn) ? false : true;

//printf("\nREAD LPA: %u\n", lpa);

    if(rb_or_gc == true) { // Sequential lpa binding
        int lpa_locate = lpa % PAGE_PER_BLOCK;

        // Get sequential possible delta first
        if(lpa_locate != 0) {
            for(int i=lpa_locate-1; i>=0; i--) {
                if(global_symb.lpa_flag[pbn][i] == 0) {
                    break;
                }
                else {
                    max_seq_cnt++;
                }
            }
        }

//printf("max_seq_cnt: %d\n", max_seq_cnt);

        if(max_seq_cnt == 0) {
            hashkey = hashing_key(lpa);

            for(int idx=ppn-1, j=new_bf_idx; idx>0; idx--) {
                if(global_symb.ppa_flag[pbn][idx] == 1) {
                    if(symbol_check(bf_man->bits_per_pg[j], j, hashkey + j, \
                                global_symb.symbol[chip][blk], st_man->sym_start, \
                                st_man->sym_bits_pg) == true) {
                        if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == lpa) {
                            found_cnt++;
                            return;
                        }
                        else {
                            notfound_cnt++;
                        }
                    }
                    else {
                        false_cnt++;
                    }

                    j--;
                }

                read_loop++;
            }
        }
        else {
            int seq_cnt = 0;
            for(int idx=ppn-1, j=new_bf_idx; idx>0; idx--) {
                if(global_symb.ppa_flag[pbn][idx] == 1) {
                    for(int delta=0; delta<=max_seq_cnt; delta++) {
                        hashkey = hashing_key(lpa - delta);

//printf("symbol checking lpa %u in new_bf_idx %d (idx %d)\n", lpa-delta, j, idx);
//printf("seq_cnt: %d\n\n", seq_cnt);

                      if(symbol_check(bf_man->bits_per_pg[j], j, hashkey + j, \
                                    global_symb.symbol[chip][blk], st_man->sym_start, \
                                    st_man->sym_bits_pg) == true) {
                            if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob ==
                                    (lpa - delta)) {
                                found_cnt++;
                                return;
                            }
                            else {
                                notfound_cnt++;
                            }
                        }
                        else {
                            false_cnt++;
                        }
                    }

                    seq_cnt = 0;
                    j--;
                }
                else {
                    seq_cnt++;
                }

                read_loop++;
            }
        }
    }
    else { // General case of bloom_read (no sequential lpa binding)
        hashkey = hashing_key(lpa);

        for(int idx=ppn-1; idx>0; idx--) {
            if(symbol_check(bf_man->bits_per_pg[idx], idx, hashkey + idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg) == true) {
                if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == lpa) {
                    found_cnt++;
                    return;
                }
                else {
                    notfound_cnt++;
                }
            }
            else {
                false_cnt++;
            }

            read_loop++;
        }
    }
    
    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[0+max_seq_cnt].oob != lpa) {
        printf("This should not happen !!\n");
        printf("LPA: %u (Until now, # of found: %d)\n", lpa, found_cnt);
        exit(1);
    }
    else {
        found_cnt++;
    }
#elif SUPERBLK
    //
#endif
}

/* 16K buffered write with Reblooming
void bloom_write(uint32_t lpa, uint32_t value, char* pttr) {
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, ppn;
    int lpa_start, buf_sz, cur_num_bf;
    Page prev_page;

    pbn = lpa >> mask; // random write
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;
    
    ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;
    
    if(ppn >= BLOCK_PPN) {
        printf("\nBefore GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tppa_flag(BF exists)\tlpa_flag(lpa exists)\tlpa_list\n");
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t\t\t%d\t\t\t%u\n", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
        }
        printf("\n");

        bloom_gc(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0);

        printf("\nAfter GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        uint64_t sumbit = 0;
        uint64_t targetsize = 0;
        for(int bf=0; bf<sblk_man->num_bf[pbn]; bf++) {
            sumbit += st_man->sym_bits_pg[bf];
        }
        sumbit += (2 * PAGE_PER_BLOCK);
        targetsize = sumbit / 8;
        if(sumbit % 8) {
            targetsize++;
        }
        printf("Sum of rebloomed ST bits in one block: [%lu bits, %lu bytes]", sumbit, targetsize);
        printf(" (%.2lf%% of PFTL)\n", (double)sumbit/PFTL*100);

        printf("ppn\tppa\twrite_lpa\tppa_flag(BF exists)\tlpa_flag(lpa exists)\tlpa_list\n");
        int sum_ppa_flag=0;
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t\t\t%d\t\t\t%u\n", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
            if(global_symb.ppa_flag[pbn][p] == 1) {
                sum_ppa_flag++;
            }
        }
        printf("\n");
        printf("GC result: # of BF: %d\n", sum_ppa_flag);
        uint64_t bfsum = 0;
        for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
            bfsum += st_man->sym_bits_pg[p];
        }
        bfsum += (2 * PAGE_PER_BLOCK);
        //uint64_t targetsize = bfsum / 8;
        targetsize = bfsum / 8;
        if(bfsum % 8) {
            targetsize++;
        }
        printf("Sum of ST bits after GC in one block: [%lu bits, %lu bytes]", bfsum, targetsize);
        printf(" (%.2lf%% of PFTL)\n", (double)bfsum/PFTL*100);
        printf("\n");
    }

    // The number of valid BFs in one superblock exceeds 50%
    if(sblk_man->num_bf[pbn] >= REBLOOM) {
        printf("\nBefore REBLOOMING in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tppa_flag(BF exists)\tlpa_flag(lpa exists)\tlpa_list\n");
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t\t\t%d\t\t\t%u\n", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
        }
        printf("\n");

        bloom_rebloom(&pbn, &chip, &way, &chnl, &blk, &ppn);
        ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

        printf("\nAfter REBLOOMING in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tppa_flag(BF exists)\tlpa_flag(lpa exists)\tlpa_list\n");
        int sum_ppa_flag=0;
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t\t\t%d\t\t\t%u\n", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
            if(global_symb.ppa_flag[pbn][p] == 1) {
                sum_ppa_flag++;
            }
        }

        printf("Reblooming result: # of BF: %d\n", sum_ppa_flag);
        uint64_t bfsum = 0;
        for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
            bfsum += st_man->sym_bits_pg[p];
        }
        bfsum += (2 * PAGE_PER_BLOCK);
        uint64_t targetsize = bfsum / 8;
        if(bfsum % 8) {
            targetsize++;
        }
        printf("Sum of ST bits after reblooming in one block: [%lu bits, %lu bytes]", bfsum, targetsize);
        printf(" (%.2lf%% of PFTL)\n", (double)bfsum/PFTL*100);
        printf("\n");
    }
    
    buf_sz = block_buffer[chip][blk].buf_sz;

    // Buffering lpa
    if(buf_sz < BUF_SZ) {
        block_buffer[chip][blk].lpa[buf_sz] = lpa;
        block_buffer[chip][blk].value[buf_sz] = value;
        block_buffer[chip][blk].buf_sz++; buf_sz++;
    }

    // Flush write buffer
    if(buf_sz == BUF_SZ) {
        lpa_start = ppn * BUF_SZ;
        cur_num_bf = sblk_man->num_bf[pbn];

        if(lpa_start != 0) {
            prev_page = storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start - 1];
        }
        
        // TODO: how to synchronize BF's start? 0 or 1?
        for(int i=0; i<BUF_SZ; i++) {
            if(lpa_start + i != 0) {
                if(prev_page.oob + 1 == block_buffer[chip][blk].lpa[i]) {
                    global_symb.ppa_flag[pbn][lpa_start + i] = 0;
                }
                else {
                    hashkey = hashing_key(block_buffer[chip][blk].lpa[i]);
                    symbol_set(bf_man->bits_per_pg[cur_num_bf + i], cur_num_bf + i, \
                            hashkey + cur_num_bf + i, global_symb.symbol[chip][blk], \
                            st_man->sym_start, st_man->sym_bits_pg, false);
//printf("write_lpa: %u\tnew_bf_idx: %d\n", block_buffer[chip][blk].lpa[i], cur_num_bf+i);

                    sblk_man->num_bf[pbn]++;
                    global_symb.ppa_flag[pbn][lpa_start + i] = 1;
                }
            }

            storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start + i].page = block_buffer[chip][blk].value[i];
            storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start + i].oob = block_buffer[chip][blk].lpa[i];
            
            global_symb.lpa_flag[pbn][block_buffer[chip][blk].lpa[i] % PAGE_PER_BLOCK] = 1;
            //global_symb.ppa_flag[pbn][lpa_start + i] = 1;
prev_page.oob = block_buffer[chip][blk].lpa[i];
            disk_write_cnt++;
        }

        storage.chip_arr[way][chnl].data_blk[blk].empty++;
        buf_sz = block_buffer[chip][blk].buf_sz = 0;
        
        write_cnt += BUF_SZ;
    }
}
*/

/* 16K buffered write with Reblooming
void bloom_read(uint32_t lpa) {
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, ppn;
    struct timeval strt, end;

    bool rb_or_gc;
    int cur_key = lpa % PAGE_PER_BLOCK;
    uint32_t real_lpa;
    int delta=0, new_bf_idx;
    
    pbn = lpa >> mask;
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;

    ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

    new_bf_idx = sblk_man->num_bf[pbn];
    rb_or_gc = (new_bf_idx == (ppn * BUF_SZ)) ? false : true;
//printf("READ REQUEST LPA %u in num_buf %d\n", lpa, new_bf_idx);

    if(rb_or_gc == true) {
        // Get sequential delta first
        if(cur_key != 0) {
            for(int idx=cur_key-1; idx>=0; idx--) {
                if(global_symb.lpa_flag[pbn][idx] == 0) {
                    break;
                }
                else {
//printf("deltachk idx: %d\n", idx);
                    delta++;
                }
            }
        }

        real_lpa = lpa - delta;
        hashkey = hashing_key(lpa - delta);
        
        int seq_cnt = 0;
        for(int idx=ppn*BUF_SZ-1, j=new_bf_idx; idx>0; idx--) {
            if(global_symb.ppa_flag[pbn][idx] == 1) {
                if(symbol_check(bf_man->bits_per_pg[j], j, hashkey + j, \
                            global_symb.symbol[chip][blk], st_man->sym_start, \
                            st_man->sym_bits_pg) == true) {
                    if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == real_lpa) {
                        found_cnt++;
                        return;
                    }
                    else {
                        notfound_cnt++;
                    }
                }
                else {
                    false_cnt++;
                }
//printf("read_lpa: %u\tidx: %d\tnew_bf_idx: %d\n", lpa, idx, j);
//printf("delta: %d\tcur_key: %d\n\n", delta, cur_key);

                j--;
            }

            read_loop++;
        }
    }
    else {
        hashkey = hashing_key(lpa);

        for(int idx=ppn*BUF_SZ-1; idx>0; idx--) {
            if(symbol_check(bf_man->bits_per_pg[idx], idx, hashkey + idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg) == true) {
                if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == lpa) {
                    found_cnt++;
                    return;
                }
                else {
                    notfound_cnt++;
                }
            }
            else {
                false_cnt++;
            }

            read_loop++;
        }
    }
    
    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[0+delta].oob != lpa) {
        printf("This should not happen !!\n");
printf("lpa: %u (", lpa);
printf("until now found_cnt: %d)\n", found_cnt);
        exit(1);
    }
    else {
        found_cnt++;
    }
}
*/

void symbol_init() {
    // Allocate SManager
    st_man = (SManager*)malloc(sizeof(SManager));
#if BLK
    st_man->sym_bits_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    st_man->new_pidx_start = (uint64_t*)malloc(sizeof(uint64_t) * TOTAL_BLOCK);
    
    memset(st_man->sym_bits_pg, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    
    for(int i=0; i<TOTAL_BLOCK; i++) {
        st_man->new_pidx_start[i] = 0;
    }
#elif SUPERBLK
    st_man->sym_bits_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SUPERBLK);
    st_man->new_pidx_start = (uint64_t*)malloc(sizeof(uint64_t) * TOTAL_SUPERBLK);

    memset(st_man->sym_bits_pg, 0, sizeof(uint64_t) * PAGE_PER_SUPERBLK);

    for(int i=0; i<TOTAL_SUPERBLK; i++) {
        st_man->new_pidx_start[i] = 0;
    }
#endif
   
    st_man->sym_bits_total = st_man->dead_bytes_total = 0;
    st_man->sym_bits_chip = st_man->sym_bits_blk = st_man->sym_bits_super_blk = 0;

#if BLK
    // Allocate symbol delimiter
    st_man->sym_start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    memset(st_man->sym_start, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    
    // Calculate symbol bits
    uint64_t symbol_bit=0, sum=0;

    for(int p=1; p<PAGE_PER_BLOCK; p++) {
#elif SUPERBLK
    // Allocate symbol delimiter
    st_man->sym_start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SUPERBLK);
    memset(st_man->sym_start, 0, sizeof(uint64_t) * PAGE_PER_SUPERBLK);

    // Calculate symbol bits
    uint64_t symbol_bit=0, sum=0;

    for(int p=1; p<PAGE_PER_SUPERBLK; p++) {
#endif
        symbol_bit = bf_man->bits_per_pg[p];
        symbol_bit = ceil(log(symbol_bit) / log(2));
		
        if(p == 1){
			symbol_bit += 4;
		}
        
        st_man->sym_bits_pg[p] = symbol_bit;
        
        st_man->sym_start[p] = sum;
        sum += symbol_bit;
    }

#if BLK
    st_man->sym_bits_blk = sum;
    st_man->sym_bits_chip = sum * BLOCK_PER_CHIP;
#elif SUPERBLK
    st_man->sym_bits_super_blk = sum;
    st_man->sym_bits_chip = sum * (SUPERBLK_PER_CHIP);
#endif
    st_man->sym_bits_total = st_man->sym_bits_chip * CHIP;
    
    st_man->sym_bytes_total = sum / 8;
    if(sum % 8) {
        st_man->sym_bytes_total++;
    }

#if BLK
    st_man->dead_bytes_total = PAGE_PER_BLOCK / 8;
#elif SUPERBLK
    st_man->dead_bytes_total = PAGE_PER_SUPERBLK / 8;
#endif

    // Allocate symbol table
    global_symb.symbol = (uint8_t***)malloc(sizeof(uint8_t**) * CHIP);
    memset(global_symb.symbol, 0, sizeof(uint8_t**) * CHIP);

    for(int c=0; c<CHIP; c++) {
#if BLK
        global_symb.symbol[c] = (uint8_t**)malloc(sizeof(uint8_t*) * BLOCK_PER_CHIP);
        memset(global_symb.symbol[c], 0, sizeof(char*) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_symb.symbol[c][b] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
            memset(global_symb.symbol[c][b], 0, sizeof(uint8_t) * st_man->sym_bytes_total);
        }
#elif SUPERBLK
        global_symb.symbol[c] = (uint8_t**)malloc(sizeof(uint8_t*) * SUPERBLK_PER_CHIP);
        memset(global_symb.symbol[c], 0, sizeof(char*) * SUPERBLK_PER_CHIP);

        for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
            global_symb.symbol[c][sb] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
            memset(global_symb.symbol[c][sb], 0, sizeof(uint8_t) * st_man->sym_bytes_total);
        }
#endif
    }

#if BLK
    global_symb.lpa_flag = (int**)malloc(sizeof(int*) * TOTAL_BLOCK);
    global_symb.ppa_flag = (int**)malloc(sizeof(int*) * TOTAL_BLOCK);
    for(int b=0; b<TOTAL_BLOCK; b++) {
        global_symb.lpa_flag[b] = (int*)malloc(sizeof(int) * PAGE_PER_BLOCK);
        global_symb.ppa_flag[b] = (int*)malloc(sizeof(int) * PAGE_PER_BLOCK);

        memset(global_symb.lpa_flag[b], 0, sizeof(int) * PAGE_PER_BLOCK);
        memset(global_symb.ppa_flag[b], 0, sizeof(int) * PAGE_PER_BLOCK);
    }
#elif SUPERBLK
    global_symb.lpa_flag = (int**)malloc(sizeof(int*) * TOTAL_SUPERBLK);
    global_symb.ppa_flag = (int**)malloc(sizeof(int*) * TOTAL_SUPERBLK);
    for(int sb=0; sb<TOTAL_SUPERBLK; sb++) {
        global_symb.lpa_flag[sb] = (int*)malloc(sizeof(int) * PAGE_PER_SUPERBLK);
        global_symb.ppa_flag[sb] = (int*)malloc(sizeof(int) * PAGE_PER_SUPERBLK);

        memset(global_symb.lpa_flag[sb], 0, sizeof(int) * PAGE_PER_SUPERBLK);
        memset(global_symb.ppa_flag[sb], 0, sizeof(int) * PAGE_PER_SUPERBLK);
    }
#endif
}

void symbol_destroy() {
    free(st_man->sym_bits_pg);
    free(st_man->new_pidx_start);
    free(st_man->sym_start);
    free(st_man);

    for(int c=0; c<CHIP; c++) {
#if BLK
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            free(global_symb.symbol[c][b]);
        }
#elif SUPERBLK
        for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
            free(global_symb.symbol[c][sb]);
        }
#endif

        free(global_symb.symbol[c]);
    }

#if BLK
    for(int b=0; b<TOTAL_BLOCK; b++) {
        free(global_symb.lpa_flag[b]);
        free(global_symb.ppa_flag[b]);
    }
#elif SUPERBLK
    for(int sb=0; sb<TOTAL_SUPERBLK; sb++) {
        free(global_symb.lpa_flag[sb]);
        free(global_symb.ppa_flag[sb]);
    }
#endif

    free(global_symb.symbol);
    free(global_symb.lpa_flag);
    free(global_symb.ppa_flag);
}

void print_stats(char* w_type, char* r_type) {
    uint64_t sum = 0;
    uint64_t targetsize = 0;
    
    printf("\nTEST DATASET: %d\n", DATA_SET);
    printf("TEST TYPE: %s write, %s read\n", w_type, r_type);

    printf("\n### BENCHMARK RESULTS ###\n");
    printf("Total read requests: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("Total false num: %d\n", false_cnt);
    printf("RAF: %.2f\n\n", (double)(found_cnt + notfound_cnt) / read_cnt);
   
    printf("Total write requests: %d\n", write_cnt);
    printf("Total number of disk write: %d\n", disk_write_cnt);
    printf("WAF: %.2f\n", (double)disk_write_cnt / write_cnt);

    printf("\n### BLOOMFILTER INFO ###\n");
    printf("Sum of BF assigned bits in all pages in 1 block: ");
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        sum += bf_man->bits_per_pg[p];
    }
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("[%lu bits, %lu bytes] (%.2lf%% of PFTL)\n", sum, targetsize, (double)sum/PFTL*100);

    uint64_t bits_per_chip = 0;
    printf("Sum of BF assigned bits in all blocks in 1 chip: ");
    bits_per_chip += (sum * BLOCK_PER_CHIP);
    targetsize = bits_per_chip / 8;
    if(bits_per_chip % 8) {
        targetsize++;
    }
    printf("[%lu bits, %lu bytes]\n", bits_per_chip, targetsize);

    uint64_t total_bits = 0;
    printf("All of BF assigned: ");
    total_bits += (bits_per_chip * CHIP);
    targetsize = total_bits / 8;
    if(total_bits % 8) {
        targetsize++;
    }
    printf("[%lu bits, %lu bytes]\n", total_bits, targetsize);
    
    printf("\n### SYMBOL TABLE INFO ###\n");
    targetsize = st_man->sym_bits_blk / 8;
    if(st_man->sym_bits_blk % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all pages in 1 block: [%lu bits, %lu bytes]", st_man->sym_bits_blk, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)st_man->sym_bits_blk/PFTL*100);

    targetsize = st_man->sym_bits_chip / 8;
    if(st_man->sym_bits_chip % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all blocks in 1 chip: [%lu bits, %lu bytes]\n", st_man->sym_bits_chip, targetsize);

    sum = st_man->sym_bits_chip * CHIP;
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Sum of ST assigned bits in all chips: [%lu bits, %lu bytes]\n", sum, targetsize);

    targetsize = st_man->sym_bits_total / 8;
    if(st_man->sym_bits_total % 8) {
        targetsize++;
    }
    printf("ALL of ST assigned: %lu bits, %lu bytes\n", st_man->sym_bits_total, targetsize);

    sum = 0;

    printf("\n### AFTER REBLOOMING ###\n");
    for(int p=0; p<sblk_man->num_bf[0]; p++) {
        sum += st_man->sym_bits_pg[p];
    }
    sum += (2 * PAGE_PER_BLOCK);
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Sum of rebloomed ST bits in one block: [%lu bits, %lu bytes]", sum, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)sum/PFTL*100);

    // TODO: fix it later !
    uint64_t sym_bits_per_chip = 0;
    printf("Sum of Rebloomed bits in all blocks in 1 chip: ");
    for(int b=0; b<BLOCK_PER_CHIP; b++) {
        for(int p=0; p<sblk_man->num_bf[b]; p++) {
            sym_bits_per_chip += st_man->sym_bits_pg[p];
        }
        sym_bits_per_chip += (2 * PAGE_PER_BLOCK);
    }
    targetsize = sym_bits_per_chip / 8;
    if(sym_bits_per_chip % 8) {
        targetsize++;
    }
    printf("[%lu bits, %lu bytes]\n", sym_bits_per_chip, targetsize);

    uint64_t total_sym_bits = 0;
    printf("Amount of SYMOL_BITS assigned: ");
    for(int w=0; w<WAY; w++) {
        for(int chnl=0; chnl<CHANNEL; chnl++) {
            int pbn = w * CHANNEL + chnl;

            for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
                total_sym_bits += st_man->sym_bits_pg[p];
            }
            total_sym_bits += (2 * PAGE_PER_BLOCK);
        }
    }
    targetsize = total_sym_bits / 8;
    if(total_sym_bits % 8) {
        targetsize++;
    }
    printf("[%lu bits, %lu bytes]\n", total_sym_bits, targetsize);

    // Time records
	printf("\nTotal write time: %.f (us)\n", w_time);
	printf("Total read time: %.f (us)\n", r_time);
	printf("Average write time: %f (us)\n", w_time/DATA_SET);
	printf("Average read time: %f (us)\n", r_time/DATA_SET);

	printf("Total read loop count: %d\n", read_loop);
	printf("Total read check time: %.f (us)\n", read_check);
	printf("Read check time per req: %f (us)\n", read_check/DATA_SET);
	printf("Average read check time: %f (us)\n", read_check/read_loop);

    printf("\nTEST COMPLETE !!\n");
    fflush(stdout);
}

int main(int argc, char** argv) {
    int try_rd = DATA_SET, try_wr = DATA_SET;
    uint32_t lpa, pbn;
    double ptv_p=0.0, false_p=0.0;
    uint32_t *read_arr, *write_arr; // LBA array
    uint8_t op;
    struct timeval strt, end;

    read_cnt = write_cnt = true_cnt = false_cnt = found_cnt = notfound_cnt = 0;

    if(argc != 3) {
        printf("Usage: ./simulationFTL W-type R-type\n(Type must be either SEQ or RAND)\n");
        exit(1);
    }

    bloom_init();
    symbol_init();

    // Allocate LBA array
    read_arr = (uint32_t*)malloc(sizeof(uint32_t) * DATA_SET);
    write_arr = (uint32_t*)malloc(sizeof(uint32_t) * DATA_SET);

    printf("\nTEST START !!\n");
    srand((unsigned int)time(NULL));

    // Generate lpa list
    // Options: W_UNIQUE | R_UNIQUE | R_RD_ORDER | R_RD_GEN
    op = R_RD_GEN;
    //op = W_UNIQUE | R_UNIQUE;
    make_test_set(write_arr, read_arr, argv[1], argv[2], op, DATA_SET);

    val = 0;
    while(try_wr--) {
        lpa = write_arr[val];

        gettimeofday(&strt, NULL);
        bloom_write(lpa, val, argv[1]);
        gettimeofday(&end, NULL);
        w_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));

        val++;
    }

/**/
printf("\nWRITE DONE\n");
#if BLK
for(int way=0; way<WAY; way++) {
    for(int chnl=0; chnl<CHANNEL; chnl++) {
        for(int blk=0; blk<BLOCK_PER_CHIP; blk++) {
            int sum_ppa_flag=0;
            int pbn = (way*CHANNEL+chnl)*BLOCK_PER_CHIP;
            printf("[Block %d]\n", pbn);
            printf("idx\twrite_lpa\tppa_flag(BF)\tlpa_flag\tlpa\t\tbf_idx\n");

            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                printf("%d\t\t%u\t\t\t", p, \
                        storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
                printf("%d\t\t\t\t%u\t\t%u\t\t\t", global_symb.ppa_flag[pbn][p], global_symb.lpa_flag[pbn][p], pbn*PAGE_PER_BLOCK+p);
                
                if(global_symb.ppa_flag[pbn][p] == 1) {
                    sum_ppa_flag++;
                    printf("%d\n", sum_ppa_flag);
                } else {
                    printf("\n");
                }
            } printf("\n");
            printf("block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        }
    }
}
#elif SUPERBLK
for(int w=0; w<WAY; w++) {
    for(int ch=0; ch<CHANNEL; ch++) {
        for(int sb=0; sb<SUPERBLK_PER_CHIP; sb++) {
            printf("[SUPERBLOCK %d in chip %d]\n", sb, w*CHANNEL+ch);
            int sum_num_bf = 0;
            int superblk = (w*CHANNEL+ch)*SUPERBLK_PER_CHIP + sb;

            for(int b=0; b<SUPERBLK_SZ; b++) {
                printf("[Block %d]\n", b);
                printf("idx\twrite_lpa\tppa_flag(BF)\tlpa_flag\tlpalist\t\tnew_bf_idx\n");

                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    printf("%d\t\t%u\t\t\t", p, storage.chip_arr[w][ch].data_blk[sb][b].page_arr[p].oob);
                    printf("%d\t\t\t\t%u\t\t%u\t\t\t", global_symb.ppa_flag[superblk][p], \
                            global_symb.lpa_flag[superblk][p], \
                            superblk*PAGE_PER_SUPERBLK+b*PAGE_PER_BLOCK+p);

                    if(global_symb.ppa_flag[superblk][p] == 1) {
                        sum_num_bf++;
                        printf("%d\n", sum_num_bf);
                    }
                    else {
                        printf("\n");
                    }
                } printf("\n");
            } printf("\n");

            printf("Superblock %d (# of valid BFs: %d)\n", superblk, sblk_man->num_bf[superblk]);
        } printf("\n");
    }
}
fflush(stdout);
#endif
/**/
free(read_arr); free(write_arr);
return 0;

    val = 0;
    while(try_rd--) {
        lpa = read_arr[val];

        gettimeofday(&strt, NULL);
        bloom_read(lpa);
        gettimeofday(&end, NULL);
        r_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));

        val++; read_cnt++;
    }
    
    print_stats(argv[1], argv[2]);

    free(read_arr);
    free(write_arr);

    symbol_destroy();
    bloom_destroy();

    return 0;
}
