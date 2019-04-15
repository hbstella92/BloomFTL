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

// Symbol buffer
#define SYMBOL_CHUNK 10

// Options for lba generation
#define W_UNIQUE 0x1
#define R_UNIQUE 0x2
#define R_RD_ORDER 0x4
#define R_RD_GEN 0x8

// Memory prefetch
#ifdef _PREF
#define CL 64
#define TL 2
#define ML 211
#endif

double read_check=0.0;
int read_loop=0;
double w_time=0.0;
double r_time=0.0;

int data_sz;
uint32_t val;

uint32_t mask;

typedef struct {
    uint32_t lbanum;
    int vbit;
    int level;
    int found;
    int notfound;
    double fp;
} Dat;

// Struct 'Dat' follows generated sequences of their own
Dat written[DATA_SET];
Dat reading[DATA_SET]; // For statistics

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

// Symbol-related
typedef struct {
    uint8_t*** symbol;
} GlobalSymb;

GlobalSymb global_symb;

typedef struct {
    uint64_t sym_bits_chip; // 1 chip
    uint64_t sym_bits_blk; // 1 block
    uint64_t* sym_bits_pg;
    uint64_t sym_bits_total;
    uint64_t* start; // Unit: bit
} SManager;

// Symbol-related
SManager* st_man;

void bloom_init() {
    double true_p=0.0, false_p=0.0;

    mask = (1 << (int)(log(TOTAL_BLOCK) / log(2))) - 1;
   
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
                storage.chip_arr[w][ch].data_blk[b].page = (uint32_t*)malloc(sizeof(uint32_t) * PAGE_PER_BLOCK);
                storage.chip_arr[w][ch].data_blk[b].oob = (uint32_t*)malloc(sizeof(uint32_t) * PAGE_PER_BLOCK);
        
                // Initialize all values in block
                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    storage.chip_arr[w][ch].data_blk[b].page[p] = storage.chip_arr[w][ch].data_blk[b].oob[p] = 99999;
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
            
            bf_free(global_bf[c][b].bfchip_arr, PAGE_PER_BLOCK);
        }

        free(bf_man->bytes_arr[c]);
        free(bf_man->bits_per_pg[c]);
        free(bf_man->bits_per_blk[c]);

        free(global_bf[c]);
    }

    free(bf_man->bytes_arr);
    free(bf_man->bits_per_pg);
    free(bf_man->bits_per_blk);
    free(bf_man->bits_per_chip);
    free(bf_man);
    free(global_bf);

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

// W_UNIQUE works in RAND
// R_UNIQUE works in RAND (If W_UNIQUE is not set, R_UNIQUE doesn't work)
// R_RD_ORDER works in RAND READ	
// R_RD_GEN works in RAND READ.
void make_test_set(uint32_t *w_arr, uint32_t *r_arr, char *w_t, char *r_t, uint8_t options, uint32_t test_size) {
	uint8_t *page_usage;
	uint32_t make_cnt, lba, pbn;
	uint8_t read = 1, write = 1;

	if(!strcmp(w_t, "RAND")) {
		write = 0;
	}
	if(!strcmp(r_t, "RAND")) {
		read = 0;
	}

	// Generate write lba
	if(!write && !(options & W_UNIQUE)) {
		// Check whether the physical block corresponding to lba is full or not (Because there is no GC algorithm)
		make_cnt = 0;
		page_usage = (uint8_t*)malloc(TOTAL_PAGE);
		memset(page_usage, 0, TOTAL_PAGE);

		while(make_cnt < test_size) {
			lba = rand() % TOTAL_PAGE;
			pbn = lba & mask;

			if(page_usage[pbn] == PAGE_PER_BLOCK) {
				continue;
			}

			w_arr[make_cnt++] = lba;
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

	// Generate read lba
	if(read) { // Seq read
		memcpy(r_arr, w_arr, sizeof(uint32_t) * test_size);
		if(!write) {
			qsort(r_arr, test_size, sizeof(uint32_t), compare);
		}
	}
	else { // Rand read
		if((write && !(options & R_UNIQUE)) || (!write && ((options & (W_UNIQUE | R_UNIQUE)) == W_UNIQUE ||\
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

/*
uint32_t generate_lba(char req, char* pttr) {
    if(req == 'W') {        
        if(!strcmp(pttr, "SEQ")) {
            return val;
        }
        else {
            // Random and duplicated
            return rand() % DATA_SET;

            // Random and uniform
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
*/

// For debugging
void print_byte_as_bit(unsigned char* val, size_t bytes) {
    for(int i=0; i<bytes; i++) {
        for(int j=7; j>=0; j--) {
            printf("%c", (val[i] & (1 << j)) ? '1' : '0');
        } printf(" ");
    } printf("\n");
}

void bloom_write(uint32_t lba, uint32_t value, char* pttr) {
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, offset;

    pbn = (lba & mask);
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;
    offset = storage.chip_arr[way][chnl].data_blk[blk].empty;
   
    if(offset >= PAGE_PER_BLOCK) {
        printf("Generated write lba is out of range\n");
        exit(1);
    }
    
    // Do not have to make bfchip_arr of page 0
    if(offset != 0) {
        hashkey = hashing_key(lba);
		symbol_set(global_bf[chip][blk].bfchip_arr, offset, hashkey + offset, \
				global_symb.symbol[chip][blk], st_man->start, st_man->sym_bits_pg);
    }

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
    uint32_t pbn, hashkey;
    int chip, way, chnl, blk, offset;
#ifdef _PREF
    char *next_p, *end_p;
#endif
    struct timeval strt, end;    
    
    pbn = (lba & mask);
    chip = pbn / BLOCK_PER_CHIP;
    way = chip / CHANNEL;
    chnl = chip % CHANNEL;
    blk = pbn % BLOCK_PER_CHIP;
    offset = storage.chip_arr[way][chnl].data_blk[blk].empty;
#ifdef _PREF
    next_p = global_bf[chip][blk].bfchip_arr[0]->body + CL;
    end_p = global_bf[chip][blk].bfchip_arr[0]->body + ML;
#endif
    
    hashkey = hashing_key(lba);

    for(int idx=offset-1; idx>0; idx--) {
        // Check BF (Modified with Symbol table)
		if(symbol_check(global_bf[chip][blk].bfchip_arr, idx, hashkey + idx, \
					global_symb.symbol[chip][blk], st_man->start, st_man->sym_bits_pg) == true) {
            if(storage.chip_arr[way][chnl].data_blk[blk].oob[idx] == lba) {
                found_cnt++;
                reading[read_cnt].lbanum = lba;
                reading[read_cnt].level = offset - idx;
                reading[read_cnt].found++;
                return;
            }
            else {
                notfound_cnt++;
                reading[read_cnt].notfound++;
                continue;
            }
        }
        else {
            false_cnt++;
        }

#ifdef _PREF
        for(; next_p<end_p; next_p+=CL) {
            __builtin_prefetch(next_p, 0, TL);
        }
#endif
        
        read_loop++;
    }

    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[blk].oob[0] != lba) {
        printf("This should not happen !!\n");
        exit(1);
    }
    else {
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
    st_man->start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_BLOCK);
    memset(st_man->start, 0, sizeof(uint64_t) * PAGE_PER_BLOCK);
    
    // Calculate symbol bits
    uint64_t symbol_bit=0, sum=0;

    for(int p=1; p<PAGE_PER_BLOCK; p++) {
        symbol_bit = bf_man->bits_per_pg[0][0][p];
        symbol_bit = ceil(log(symbol_bit) / log(2));
		if(p == 1){
			symbol_bit += 4;
		}
        st_man->sym_bits_pg[p] = symbol_bit;
        
        st_man->start[p] = sum;
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
    global_symb.symbol = (uint8_t***)malloc(sizeof(uint8_t**) * CHIP);
    memset(global_symb.symbol, 0, sizeof(uint8_t**) * CHIP);

    for(int c=0; c<CHIP; c++) {
        global_symb.symbol[c] = (uint8_t**)malloc(sizeof(uint8_t*) * BLOCK_PER_CHIP);
        memset(global_symb.symbol[c], 0, sizeof(char*) * BLOCK_PER_CHIP);

        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_symb.symbol[c][b] = (uint8_t*)malloc(sizeof(uint8_t) * symbol_sum_byte);
            memset(global_symb.symbol[c][b], 0, sizeof(uint8_t) * symbol_sum_byte);
        }
    }
}

void symbol_destroy() {
    free(st_man->sym_bits_pg);
    free(st_man->start);
    free(st_man);

    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            free(global_symb.symbol[c][b]);
        }

        free(global_symb.symbol[c]);
    }

    free(global_symb.symbol);
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

#ifdef _PREF
    printf("\n### Memory prefetch for read ###\n");
    switch(TL) {
        case 0:
            printf("Temporal locality - none: 0\n");
            break;
        case 1:
            printf("Temporal locality - low: 1\n");
            break;
        case 2:
            printf("Temporal locality - medium: 2\n");
            break;
        case 3:
            printf("Temporal locality - high: 3\n");
            break;
    }
#endif

    // Time records
	printf("Total write time: %.f (us)\n", w_time);
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
    uint32_t lba, pbn;
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

    printf("TEST START !!\n");
    srand((unsigned int)time(NULL));

    // Generate lba list
    // Options: W_UNIQUE | R_UNIQUE | R_RD_ORDER | R_RD_GEN
    op = R_RD_ORDER;
    make_test_set(write_arr, read_arr, argv[1], argv[2], op, DATA_SET);
    //make_test_set(write_arr, read_arr, write_t, read_t, op, DATA_SET);

    val = data_sz = 0;
    while(try_wr--) {
        lba = write_arr[val];

        gettimeofday(&strt, NULL);
        bloom_write(lba, val, argv[1]);
        gettimeofday(&end, NULL);
        w_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));
        
        val++; write_cnt++;
    }

    val = 0;
    while(try_rd--) {
        lba = read_arr[val];

        gettimeofday(&strt, NULL);
        bloom_read(lba);
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
