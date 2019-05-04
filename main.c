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

#define CHANNEL 2
#define WAY 2
#define CHIP ((CHANNEL)*(WAY))

#define BLOCK_PER_CHIP 1
#define PAGE_PER_BLOCK 512 // logical unit
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP))
#define TOTAL_PAGE ((TOTAL_BLOCK)*(PAGE_PER_BLOCK))
#define LP_SZ (4*(K)) // logical page size
#define CAPACITY ((LP_SZ)*(TOTAL_PAGE))
#define PFTL ((PAGE_PER_BLOCK)*32)

//#define DATA_SET (TOTAL_PAGE*2)
#define DATA_SET (TOTAL_PAGE)

// Buffered write
#define BUF_SZ 4
#define PP_SZ ((LP_SZ)*(BUF_SZ)) // physical page size
#define BLOCK_PPN ((PAGE_PER_BLOCK)/(BUF_SZ)) // 128

// Reblooming-related
#define REBLOOM (int)((PAGE_PER_BLOCK)*0.5)

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

typedef struct {
    uint32_t lpanum;
    int vbit;
    int level;
    double fp;
} Dat;

// Struct 'Dat' follows generated sequences of their own
Dat written[DATA_SET];
Dat reading[DATA_SET]; // for statistics

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

typedef struct {
    uint32_t* lpa;
    uint32_t* value;
    int buf_sz;
} BlockBuffer; // TODO: change to struct Page

BlockBuffer** block_buffer;

typedef struct {
    uint32_t page;
    uint32_t oob;
} Page;

typedef struct {
    Page* page_arr;
    int empty; // offset of physical unit
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
    uint64_t* sym_bits_pg;
    uint64_t sym_bits_total;
    uint64_t* new_pidx_start;

    uint64_t* sym_start; // bit unit
} SManager;

SManager* st_man;

void bloom_init() {
    double true_p=0.0, false_p=0.0;

    mask = (int)(log(PAGE_PER_BLOCK) / log(2));
    //mask = (1 << (int)(log(PAGE_PER_BLOCK) / log(2))) - 1;

    // Alloc && Initialize SBlkManager
    sblk_man = (SBlkManager*)malloc(sizeof(SBlkManager));
    sblk_man->num_bf = (int*)malloc(sizeof(int) * TOTAL_BLOCK);
    for(int i=0; i<TOTAL_BLOCK; i++) {
        sblk_man->num_bf[i] = 0;
    }
   
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

    // Alloc && Initialize global_bf
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

    // Alloc memory for all blocks and pages
    storage.chip_arr = (Chip**)malloc(sizeof(Chip*) * WAY);
    for(int w=0; w<WAY; w++) {
        storage.chip_arr[w] = (Chip*)malloc(sizeof(Chip) * CHANNEL);

        for(int ch=0; ch<CHANNEL; ch++) {
            storage.chip_arr[w][ch].data_blk = (Block*)malloc(sizeof(Block) * BLOCK_PER_CHIP);

            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                storage.chip_arr[w][ch].data_blk[b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
     
                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    storage.chip_arr[w][ch].data_blk[b].page_arr[p].page = 
                        storage.chip_arr[w][ch].data_blk[b].page_arr[p].oob = 99999;
                }

                storage.chip_arr[w][ch].data_blk[b].empty = 0;
            }
        }
    }

    // Alloc memory for block buffer
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
}

void bloom_destroy() {
    free(sblk_man->num_bf);
    free(sblk_man);

    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            free(bf_man->bytes_arr[c][b]);
            free(bf_man->bits_per_pg[c][b]);
            
            bf_free(global_bf[c][b].bfchip_arr, PAGE_PER_BLOCK);

            free(block_buffer[c][b].lpa);
            free(block_buffer[c][b].value);
        }

        free(bf_man->bytes_arr[c]);
        free(bf_man->bits_per_pg[c]);
        free(bf_man->bits_per_blk[c]);

        free(global_bf[c]);

        free(block_buffer[c]);
    }

    free(bf_man->bytes_arr);
    free(bf_man->bits_per_pg);
    free(bf_man->bits_per_blk);
    free(bf_man->bits_per_chip);
    free(bf_man);
    free(global_bf);

    free(block_buffer);

    for(int w=0; w<WAY; w++) {
        for(int ch=0; ch<CHANNEL; ch++) {
            for(int b=0; b<BLOCK_PER_CHIP; b++) {
                free(storage.chip_arr[w][ch].data_blk[b].page_arr);
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
            pbn = lpa >> mask;
            //pbn = lpa & mask; // original!

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

int lpa_compare(void* a, void* b) {
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

void symbol_resymbolize(uint32_t pbn, int chip, int way, int chnl, int blk, int valid_start, int num_flush) {
    uint32_t hashkey;

    storage.chip_arr[way][chnl].data_blk[blk].empty += (num_flush / BUF_SZ);
    
    free(global_symb.symbol[chip][blk]);
    sblk_man->num_bf[pbn] = 0;

    global_symb.symbol[chip][blk] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
    memset(global_symb.symbol[chip][blk], 0, sizeof(uint8_t) * st_man->sym_bytes_total);

    int ppa_end_idx = storage.chip_arr[way][chnl].data_blk[blk].empty * BUF_SZ;
    int new_bf_idx = 1;

    // Set new symbol and Set dead bit for lpa and ppa (TODO: BF's new start idx is 1)
    for(int i=0; i<ppa_end_idx; i++) {
        if(i < valid_start) {
            if(global_symb.ppa_flag[pbn][i] == 1) {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[chip][blk][new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);

                global_symb.lpa_flag[pbn][i] = 1;

                new_bf_idx++;
            }
        }
        else if(i > valid_start) {
            if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i-1].oob + 1 == 
                    storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob) {
                global_symb.lpa_flag[pbn][i] = 0;
            }
            else {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[chip][blk][new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);

                global_symb.lpa_flag[pbn][i] = 1;

                new_bf_idx++;
            }

            global_symb.ppa_flag[pbn][i] = 1;
        }
        else {
            if(i == 0) {
                hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                symbol_set(bf_man->bits_per_pg[chip][blk][new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);

                global_symb.lpa_flag[pbn][i] = 1;

                new_bf_idx++;
            }
            else {
                if(global_symb.ppa_flag[pbn][i-1] == 1) {
                    if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i-1].oob + 1 == 
                            storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob) {
                        global_symb.lpa_flag[pbn][i] = 0;
                    }
                    else {
                        hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                        symbol_set(bf_man->bits_per_pg[chip][blk][new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                                global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);

                        global_symb.lpa_flag[pbn][i] = 1;

                        new_bf_idx++;
                    }
                }
                else {
                    hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[blk].page_arr[i].oob);
                    symbol_set(bf_man->bits_per_pg[chip][blk][new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
                            global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);

                    global_symb.lpa_flag[pbn][i] = 1;

                    new_bf_idx++;
                }
            }
            
            global_symb.ppa_flag[pbn][i] = 1;
        }
    }

    sblk_man->num_bf[pbn] = new_bf_idx - 1;
}

void bloom_gc(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, Page* evict_list, int num_evict) {
    Page *candi_list, *valid_list, prev_page;
    uint32_t hashkey;
    int num_valid = 0, candi_sz = 0;

    candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
    valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

    // Copy valid pages in evict block
    for(int i=0; i<PAGE_PER_BLOCK; i++) {
        if(global_symb.ppa_flag[*pbn][i] == 1) {
            candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[i];
        }
    }

    if(evict_list != NULL) {
        memcpy(&candi_list[candi_sz], evict_list, sizeof(Page) * num_evict);
        candi_sz += num_evict;
    }

    // Sort valid pages with lpa order
    qsort(candi_list, candi_sz, sizeof(Page), lpa_compare);

/*
printf("GC candi_sz: %d\n", candi_sz);
for(int i=0; i<candi_sz; i++) {
    printf("idx: %d\tlpa_in_candi_lst: %u\n", i, candi_list[i].oob);
} printf("\n\n");
*/

    // Free exist block
    free(storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr);

    // Re-allocate and Initialize new block
    storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p].page = 
            storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p].oob = 99999;
    }
    storage.chip_arr[*way][*chnl].data_blk[*blk].empty = 0;
    
    memset(global_symb.lpa_flag[*pbn], -1, sizeof(int) * PAGE_PER_BLOCK);
    memset(global_symb.ppa_flag[*pbn], -1, sizeof(int) * PAGE_PER_BLOCK);

    // Elect single valid pages from candi_list
    prev_page = valid_list[0] = candi_list[0];
    for(int i=1; i<candi_sz; i++) {
        if(prev_page.oob != candi_list[i].oob) {
            num_valid++;
        }

        valid_list[num_valid] = candi_list[i];
        prev_page = candi_list[i];
    }
    num_valid++;

/*
printf("GC num_valid: %d\n", num_valid);
for(int i=0; i<num_valid; i++) {
    printf("idx: %d\tlpa_in_valid_lst: %u\n", i, valid_list[i].oob);
} printf("\n\n");
*/

    int num_need_buf = num_valid % BUF_SZ;
    bool need_buf = num_need_buf ? true : false;

    // Flush valid pages
    memcpy(storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr, valid_list, sizeof(Page) * (num_valid - num_need_buf));
    symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, 0, (num_valid - num_need_buf));
    
    *ppn = storage.chip_arr[*way][*chnl].data_blk[*blk].empty;
    disk_write_cnt += (num_valid - num_need_buf);

    // Buffer if there are remained valid pages
    if(need_buf == true) {
        int buf_sz = block_buffer[*chip][*blk].buf_sz;

        for(int i=buf_sz, j=0; i<buf_sz+num_need_buf; i++, j++) {
            if(i == BUF_SZ) {
                int lpa_start = (*ppn) * BUF_SZ;
                
                if(lpa_start == PAGE_PER_BLOCK) {
                    printf("This is GC error\n");
                    exit(1);
                }

                for(int k=0; k<BUF_SZ; k++) {
                    // TODO: how to synchronize BF's start?
                    if(lpa_start + k != 0) {
                        hashkey = hashing_key(block_buffer[*chip][*blk].lpa[k]);
                        symbol_set(bf_man->bits_per_pg[*chip][*blk][lpa_start + k], lpa_start + k, hashkey + lpa_start + k, \
                                global_symb.symbol[*chip][*blk], st_man->sym_start, st_man->sym_bits_pg, false);

                        sblk_man->num_bf[*pbn]++;
                    }

                    storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].page = block_buffer[*chip][*blk].value[k];
                    storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].oob = block_buffer[*chip][*blk].lpa[k];

                    global_symb.ppa_flag[*pbn][lpa_start + k] = 1;
                    global_symb.lpa_flag[*pbn][lpa_start + k] = 1;

                    disk_write_cnt++;
                }

                storage.chip_arr[*way][*chnl].data_blk[*blk].empty++;
                buf_sz = block_buffer[*chip][*blk].buf_sz = 0;
            }

            block_buffer[*chip][*blk].lpa[(i % BUF_SZ)] = valid_list[num_valid - num_need_buf + j].oob;
            block_buffer[*chip][*blk].value[(i % BUF_SZ)] = valid_list[num_valid - num_need_buf + j].page;
            block_buffer[*chip][*blk].buf_sz++;
        }

/*
printf("BUFFER UNFINISHED DATA\n");
for(int i=0; i<block_buffer[*chip][*blk].buf_sz; i++) {
    printf("buffer_idx: %d\tbuffer_lpa: %u\n", i, block_buffer[*chip][*blk].lpa[i]);
}
*/

    }
    
    free(valid_list);
    free(candi_list);
}

void bloom_rebloom(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn) {
    Page *candi_list, *evict_list, prev_page;
    uint32_t hashkey;
    int num_candi = 0, num_evict = 0, candi_sz = 0;

    num_candi = storage.chip_arr[*way][*chnl].data_blk[*blk].empty * BUF_SZ;

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
                global_symb.lpa_flag[*pbn][j] = 0;
            }
        }
    }

    if(num_candi + num_evict >= PAGE_PER_BLOCK) {
        bloom_gc(pbn, chip, way, chnl, blk, ppn, evict_list, num_evict);

        free(evict_list);
        free(candi_list);
        return;
    }

/*
printf("RB num_candi: %d\n", num_candi);
for(int i=0; i<num_candi; i++) {
    printf("idx: %d\tlpa_in_candi_lst: %u\n", i, candi_list[i].oob);
} printf("\n\n");

printf("RB num_evict: %d\n", num_evict);
for(int i=0; i<num_evict; i++) {
    printf("idx: %d\tlpa_in_evict_lst: %u\n", i, evict_list[i].oob);
} printf("\n");
*/

    // Append rebloomed pages
    bool need_buf = false;
    int buf_sz;
    int remain_evict = num_evict;
    int num_need_buf = remain_evict % BUF_SZ;
    
    *ppn = storage.chip_arr[*way][*chnl].data_blk[*blk].empty;

    if(!num_need_buf) {
        for(int i=0, p=num_candi; i<num_evict; i++, p++) {
            if(p == PAGE_PER_BLOCK) {
                symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, num_evict);
                
                bloom_gc(pbn, chip, way, chnl, blk, ppn, NULL, 0);

                p = 0;
                num_candi = 0;
                remain_evict = num_evict - i;
                num_need_buf = remain_evict % BUF_SZ;

                if(num_need_buf != 0) {
                    remain_evict -= num_need_buf;
                    need_buf = true;
                }
            }

            if((p == remain_evict) && need_buf) {
                break;
            }

            storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p] = evict_list[i];

            disk_write_cnt++;
        }

        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, remain_evict);

        if(need_buf == true) {
            buf_sz = block_buffer[*chip][*blk].buf_sz;

            for(int i=buf_sz, j=0; i<buf_sz+num_need_buf; i++, j++) {
                if(i == BUF_SZ) {
                    int lpa_start = (*ppn) * BUF_SZ;

                    if(lpa_start == PAGE_PER_BLOCK) {
                        bloom_gc(pbn, chip, way, chnl, blk, ppn, NULL, 0);

                        lpa_start = (*ppn) * BUF_SZ;
                    }
                    
                    for(int k=0; k<BUF_SZ; k++) {
                    // TODO: how to synchronize BF's start?
                        if(lpa_start + k != 0) {
                            hashkey = hashing_key(block_buffer[*chip][*blk].lpa[k]);
                            symbol_set(bf_man->bits_per_pg[*chip][*blk][lpa_start + k], lpa_start + k, hashkey + lpa_start + k, \
                                    global_symb.symbol[*chip][*blk], st_man->sym_start, st_man->sym_bits_pg, false);

                            sblk_man->num_bf[*pbn]++;
                        }

                        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].page = block_buffer[*chip][*blk].value[k];
                        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].oob = block_buffer[*chip][*blk].lpa[k];

                        global_symb.ppa_flag[*pbn][lpa_start + k] = 1;
                        global_symb.lpa_flag[*pbn][lpa_start + k] = 1;

                        disk_write_cnt++;
                    }

                    storage.chip_arr[*way][*chnl].data_blk[*blk].empty++;
                    buf_sz = block_buffer[*chip][*blk].buf_sz = 0;
                }
                
                block_buffer[*chip][*blk].lpa[(i % BUF_SZ)] = evict_list[num_evict-num_need_buf+j].oob;
                block_buffer[*chip][*blk].value[(i % BUF_SZ)] = evict_list[num_evict-num_need_buf+j].page;
                block_buffer[*chip][*blk].buf_sz++;
            }

/*
printf("BUFFER UNFINISHED DATA\n");
for(int i=0; i<block_buffer[*chip][*blk].buf_sz; i++) {
    printf("buffer_idx: %d\tbuffer_lpa: %u\n", i, block_buffer[*chip][*blk].lpa[i]);
}
*/

        }
    }
    else {
        need_buf = true;

        for(int i=0, p=num_candi; i<num_evict; i++, p++) {
            if(p == PAGE_PER_BLOCK) {
                symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, num_evict);

                bloom_gc(pbn, chip, way, chnl, blk, ppn, NULL, 0);

                p = 0;
                num_candi = 0;
                remain_evict = num_evict - i;
                num_need_buf = remain_evict % BUF_SZ;

                if(num_need_buf != 0) {
                    remain_evict -= num_need_buf;
                }
                else {
                    need_buf = false;
                }
            }
            
            if(need_buf) {
                if((p == num_candi + num_evict - num_need_buf) || (p == remain_evict)) {
                    break;
                }
            }

            storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[p] = evict_list[i];

            disk_write_cnt++;
        }

        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, remain_evict);

        if(need_buf == true) {
            buf_sz = block_buffer[*chip][*blk].buf_sz;

            for(int i=buf_sz, j=0; i<buf_sz+num_need_buf; i++, j++) {
                if(i == BUF_SZ) {
                    int lpa_start = (*ppn) * BUF_SZ;

                    if(lpa_start == PAGE_PER_BLOCK) {
                        bloom_gc(pbn, chip, way, chnl, blk, ppn, NULL, 0);

                        lpa_start = (*ppn) * BUF_SZ;
                    }

                    for(int k=0; k<BUF_SZ; k++) {
                        if(lpa_start + k != 0) {
                            hashkey = hashing_key(block_buffer[*chip][*blk].lpa[k]);
                            symbol_set(bf_man->bits_per_pg[*chip][*blk][lpa_start + k], lpa_start + k, hashkey + lpa_start + k, \
                                    global_symb.symbol[*chip][*blk], st_man->sym_start, st_man->sym_bits_pg, false);

                            sblk_man->num_bf[*pbn]++;
                        }

                        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].page = block_buffer[*chip][*blk].value[k];
                        storage.chip_arr[*way][*chnl].data_blk[*blk].page_arr[lpa_start + k].oob = block_buffer[*chip][*blk].lpa[k];

                        global_symb.ppa_flag[*pbn][lpa_start + k] = 1;
                        global_symb.lpa_flag[*pbn][lpa_start + k] = 1;

                        disk_write_cnt++;
                    }

                    storage.chip_arr[*way][*chnl].data_blk[*blk].empty++;
                    buf_sz = block_buffer[*chip][*blk].buf_sz = 0;
                }
                
                block_buffer[*chip][*blk].lpa[(i % BUF_SZ)] = evict_list[num_evict-num_need_buf+j].oob;
                block_buffer[*chip][*blk].value[(i % BUF_SZ)] = evict_list[num_evict-num_need_buf+j].page;
                block_buffer[*chip][*blk].buf_sz++;
            }

/*
printf("BUFFER UNFINISHED DATA\n");
for(int i=0; i<block_buffer[*chip][*blk].buf_sz; i++) {
    printf("buffer_idx: %d\tbuffer_lpa: %u\n", i, block_buffer[*chip][*blk].lpa[i]);
}
*/

        }
    }

    free(evict_list);
    free(candi_list);
}

// 16K buffered write with Reblooming
void bloom_write(uint32_t lpa, uint32_t value, char* pttr) {
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, ppn;
    int lpa_start, buf_sz;

    pbn = lpa >> mask;
    //pbn = ((lpa >> LPA_BUF) & mask);
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;
    
    ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;
   
    if(ppn >= BLOCK_PPN) {

/*
        printf("\nGC-triggered lpa: %u\tpbn: %u\t", lpa, pbn);
        printf("chip: %d\t", chip);
        printf("way: %d\tchnl: %d\tblk: %d\n\n", way, chnl, blk); 
        printf("\nBefore GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tlpa_flag\tppa_flag(valid_bit)\n");
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t%d\n", global_symb.lpa_flag[pbn][p], global_symb.ppa_flag[pbn][p]);
        }
        printf("\n");
*/

        bloom_gc(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0);

/*
        printf("\nAfter GC in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tlpa_flag\tppa_flag(valid_bit)\n");
        int sum_lpa_flag=0;
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t%d\n", global_symb.lpa_flag[pbn][p], global_symb.ppa_flag[pbn][p]);
            if(global_symb.lpa_flag[pbn][p] == 1) {
                sum_lpa_flag++;
            }
        }
        printf("\n");
        printf("GC result: # of (lpa_flag value == 1): %d\n", sum_lpa_flag);
        uint64_t bfsum = 0;
        for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
            bfsum += st_man->sym_bits_pg[p];
        }
        uint64_t targetsize = bfsum / 8;
        if(bfsum % 8) {
            targetsize++;
        }
        printf("Sum of ST bits after GC in one block: [%lu bits, %lu bytes]", bfsum, targetsize);
        printf(" (%.2lf%% of PFTL)\n", (double)bfsum/PFTL*100);
        printf("\n");
*/

    }

    // The number of valid BFs in one superblock exceeds 50%
    if(sblk_man->num_bf[pbn] >= REBLOOM) {

/*
        printf("\nBefore REBLOOMING in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tlpa_flag\tppa_flag(valid_bit)\n");
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t%d\n", global_symb.lpa_flag[pbn][p], global_symb.ppa_flag[pbn][p]);
        }
        printf("\n");
*/

        bloom_rebloom(&pbn, &chip, &way, &chnl, &blk, &ppn);

        ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

/*        
        printf("\nAfter REBLOOMING in block %u (# of valid BFs: %d)\n", pbn, sblk_man->num_bf[pbn]);
        printf("ppn\tppa\twrite_lpa\tlpa_flag\tppa_flag(valid_bit)\n");
        int sum_lpa_flag=0;
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
            printf("%d\t\t\t%d\n", global_symb.lpa_flag[pbn][p], global_symb.ppa_flag[pbn][p]);
            if(global_symb.lpa_flag[pbn][p] == 1) {
                sum_lpa_flag++;
            }
        }
        printf("Reblooming result: # of (lpa_flag value == 1): %d\n", sum_lpa_flag);
        uint64_t bfsum = 0;
        for(int p=0; p<sblk_man->num_bf[pbn]; p++) {
            bfsum += st_man->sym_bits_pg[p];
        }
        uint64_t targetsize = bfsum / 8;
        if(bfsum % 8) {
            targetsize++;
        }
        printf("Sum of ST bits after reblooming in one block: [%lu bits, %lu bytes]", bfsum, targetsize);
        printf(" (%.2lf%% of PFTL)\n", (double)bfsum/PFTL*100);
        printf("\n");
*/

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

        // TODO: how to synchronize BF's start? 0 or 1?
        for(int i=0; i<BUF_SZ; i++) {
            if(lpa_start + i != 0) {
                hashkey = hashing_key(block_buffer[chip][blk].lpa[i]);
                symbol_set(bf_man->bits_per_pg[chip][blk][lpa_start + i], lpa_start + i, hashkey + lpa_start + i, \
                        global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg, false);
                
                sblk_man->num_bf[pbn]++;
            }

            storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start + i].page = block_buffer[chip][blk].value[i];
            storage.chip_arr[way][chnl].data_blk[blk].page_arr[lpa_start + i].oob = block_buffer[chip][blk].lpa[i];
            
            global_symb.ppa_flag[pbn][lpa_start + i] = 1;
            global_symb.lpa_flag[pbn][lpa_start + i] = 1;

            disk_write_cnt++;
        }

        storage.chip_arr[way][chnl].data_blk[blk].empty++;
        buf_sz = block_buffer[chip][blk].buf_sz = 0;
        
        write_cnt += BUF_SZ;
    }
}

// 16K buffered write without Reblooming (original backup)
void bloom_read(uint32_t lpa) {
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, offset;
    struct timeval strt, end;    
    
    pbn = lpa >> 7;
//    pbn = ((lpa >> LPA_BUF) & mask);
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;

    offset = storage.chip_arr[way][chnl].data_blk[blk].empty;
    
    hashkey = hashing_key(lpa);

    for(int idx=offset*BUF_SZ-1; idx>0; idx--) {
        if(symbol_check(bf_man->bits_per_pg[chip][blk][idx], idx, hashkey + idx, \
                    global_symb.symbol[chip][blk], st_man->sym_start, st_man->sym_bits_pg) == true) {
            if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == lpa) {
                found_cnt++;
                return;
            }
            else {
                notfound_cnt++;
                //continue;
            }
        }
        else {
            false_cnt++;
        }

        read_loop++;
    }

    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[0].oob != lpa) {
        printf("This should not happen !!\n");
        exit(1);
    }
    else {
        found_cnt++;
        reading[read_cnt].lpanum = lpa;
        reading[read_cnt].level = offset - 0;
    }
}

void symbol_init() {
    // Allocate SManager
    st_man = (SManager*)malloc(sizeof(SManager));
    st_man->sym_bits_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    st_man->new_pidx_start = (uint64_t*)malloc(sizeof(uint64_t) * TOTAL_BLOCK);
    
    memset(st_man->sym_bits_pg, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    for(int i=0; i<TOTAL_BLOCK; i++) {
        st_man->new_pidx_start[i] = 0;
    }
    st_man->sym_bits_total = st_man->dead_bytes_total = st_man->sym_bits_chip = st_man->sym_bits_blk = 0;

    // Allocate symbol delimiter
    st_man->sym_start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    memset(st_man->sym_start, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    
    // Calculate symbol bits
    uint64_t symbol_bit=0, sum=0;

    for(int p=1; p<PAGE_PER_BLOCK; p++) {
        symbol_bit = bf_man->bits_per_pg[0][0][p];
        symbol_bit = ceil(log(symbol_bit) / log(2));
		if(p == 1){
			symbol_bit += 4;
		}
        st_man->sym_bits_pg[p] = symbol_bit;
        
        st_man->sym_start[p] = sum;
        sum += symbol_bit;
    }

    st_man->sym_bits_blk = sum;
    st_man->sym_bits_chip = sum * BLOCK_PER_CHIP;
    st_man->sym_bits_total = st_man->sym_bits_chip * CHIP;
    
    st_man->sym_bytes_total = sum / 8;
    if(sum % 8) {
        st_man->sym_bytes_total++;
    }

    st_man->dead_bytes_total = PAGE_PER_BLOCK / 8;

    // Allocate symbol table
    global_symb.symbol = (uint8_t***)malloc(sizeof(uint8_t**) * CHIP);
    memset(global_symb.symbol, 0, sizeof(uint8_t**) * CHIP);

    for(int c=0; c<CHIP; c++) {
        global_symb.symbol[c] = (uint8_t**)malloc(sizeof(uint8_t*) * BLOCK_PER_CHIP);
        memset(global_symb.symbol[c], 0, sizeof(char*) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_symb.symbol[c][b] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
            memset(global_symb.symbol[c][b], 0, sizeof(uint8_t) * st_man->sym_bytes_total);
        }
    }

    global_symb.lpa_flag = (int**)malloc(sizeof(int*) * TOTAL_BLOCK);
    global_symb.ppa_flag = (int**)malloc(sizeof(int*) * TOTAL_BLOCK);
    for(int b=0; b<TOTAL_BLOCK; b++) {
        global_symb.lpa_flag[b] = (int*)malloc(sizeof(int) * PAGE_PER_BLOCK);
        global_symb.ppa_flag[b] = (int*)malloc(sizeof(int) * PAGE_PER_BLOCK);

        memset(global_symb.lpa_flag[b], -1, sizeof(int) * PAGE_PER_BLOCK);
        memset(global_symb.ppa_flag[b], -1, sizeof(int) * PAGE_PER_BLOCK);
    }
}

void symbol_destroy() {
    free(st_man->sym_bits_pg);
    free(st_man->new_pidx_start);
    free(st_man->sym_start);
    free(st_man);

    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            free(global_symb.symbol[c][b]);
        }

        free(global_symb.symbol[c]);
    }

    for(int b=0; b<TOTAL_BLOCK; b++) {
        free(global_symb.lpa_flag[b]);
        free(global_symb.ppa_flag[b]);
    }

    free(global_symb.symbol);
    free(global_symb.lpa_flag);
    free(global_symb.ppa_flag);
}

void print_stats(char* w_type, char* r_type) {
    uint64_t sum = 0;
/*
    printf("\nFOUND LEVEL OF ALL READ REQUESTS\n");
    for(int i=0; i<DATA_SET; i++) {
        printf("req_num: %d level: %d\n", i, reading[i].level);
    }
*/    
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

    printf("\n### AFTER REBLOOMING ###\n");
    sum = 0;
    for(int p=0; p<sblk_man->num_bf[0]; p++) {
        sum += st_man->sym_bits_pg[p];
    }
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Sum of rebloomed ST bits in one block: [%lu bits, %lu bytes]", sum, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)sum/PFTL*100);
/*
    sum = 0;
    for(int b=0; b<TOTAL_BLOCK; b++) {
        sum += ((st_man->sym_start[bf_man->num_bf[b] + 1]) + 128 / 8);
    }
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Sum of rebloomed ST bits in all blocks: [%lu bits, %lu bytes]\n", sum, targetsize);
    
    sum = (st_man->sym_start[bf_man->num_bf[0] + 1] + (128 / 8));
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Sum of rebloomed ST bits in one block: [%lu bits, %lu bytes]", sum, targetsize);
    printf(" (%.2lf%% of PFTL)\n", (double)sum/PFTL*100);
*/
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

/*
    val = 0;
    while(try_rd--) {
        lpa = read_arr[val];

        gettimeofday(&strt, NULL);
        bloom_read(lpa);
        gettimeofday(&end, NULL);
        r_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));

        val++; read_cnt++;
    }
*/
    print_stats(argv[1], argv[2]);

/*
printf("FINAL RESULT\n");
printf("ppn\tppa\twrite_lpa\tlpa_flag\tppa_flag(valid_bit)\n");
int sum_lpa_flag=0;
for(int way=0; way<WAY; way++) {
    for(int chnl=0; chnl<CHANNEL; chnl++) {
        for(int blk=0; blk<BLOCK_PER_CHIP; blk++) {
            printf("block %d\n", (way*CHANNEL+chnl)*BLOCK_PER_CHIP);
            
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                int pbn = (way*CHANNEL+chnl)*BLOCK_PER_CHIP;
                printf("%d\t%d\t\t%u\t\t\t", p/BUF_SZ, p, storage.chip_arr[way][chnl].data_blk[blk].page_arr[p].oob);
                printf("%d\t\t\t%d\n", global_symb.lpa_flag[pbn][p], global_symb.ppa_flag[pbn][p]);
                fflush(stdout);
                if(global_symb.lpa_flag[pbn][p] == 1) {
                    sum_lpa_flag++;
                }
            }
            printf("\n");
        }
    }
}
*/

    free(read_arr);
    free(write_arr);

    symbol_destroy();
    bloom_destroy();

    return 0;
}
