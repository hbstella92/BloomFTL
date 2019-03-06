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

//#define BLOCK_PER_CHIP 8
#define BLOCK_PER_CHIP 1
#define PAGE_PER_BLOCK 128
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP))
#define TOTAL_PAGE ((TOTAL_BLOCK)*(PAGE_PER_BLOCK))
#define PAGE_SZ (4*(K))

#define LSB (int)(log(TOTAL_BLOCK) / log(2))
#define MASK ((1 << LSB) - 1)

#define CAPACITY ((PAGE_SZ)*(TOTAL_PAGE))
#define DATA_SET (TOTAL_PAGE)

#define WRKMEM_SZ 10

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
    BF* bfchip_arr;
} GlobalBF;

typedef struct {
    uint64_t* bits_per_chip;
    uint64_t** bits_per_blk;
    uint64_t*** bits_per_pg;
    uint64_t*** bytes_arr;
} BFManager;

GlobalBF*** global_bf;
BFManager* bf_man;

Storage storage;
Block** data_block;

double* fp_arr;

// Compress && Decompress BF
typedef struct {
    char* comp_arr;
    size_t len;
} Compbuf;

Compbuf*** compbuf;
Compbuf*** decompbuf;
Compbuf*** checkbuf;

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
    global_bf = (GlobalBF***)malloc(sizeof(GlobalBF**) * CHIP);
    for(int c=0; c<CHIP; c++) {
        global_bf[c] = (GlobalBF**)malloc(sizeof(GlobalBF*) * BLOCK_PER_CHIP);
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            global_bf[c][b] = (GlobalBF*)malloc(sizeof(GlobalBF) * PAGE_PER_BLOCK);
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                if(p == 0) {
                    false_p = 0.0001;
                    global_bf[c][b][p].bfchip_arr = bf_init(1, false_p);

                    bf_man->bits_per_pg[c][b][p] = 0;
                    bf_man->bytes_arr[c][b][p] = 0;
                }
                else {
                    true_p = pow(0.9, (double)1/p);
                    false_p = 1 - true_p;
                    
                    global_bf[c][b][p].bfchip_arr = bf_init(1, false_p);

                    bf_man->bits_per_pg[c][b][p] = bf_bits(1, false_p);

                    int targetsize = bf_man->bits_per_pg[c][b][p] / 8;
                    if(bf_man->bits_per_pg[c][b][p] % 8) {
                        targetsize++;
                    }
                    bf_man->bytes_arr[c][b][p] = targetsize;
//printf("bytes: %d\n", targetsize);
                }
//                printf("entry #: %d\tfalse p: %f\tbit #: %d\tfunc #: %u\n", \
                        global_bf[c][b][p].bfchip_arr->n, \
                        global_bf[c][b][p].bfchip_arr->p, \
                        global_bf[c][b][p].bfchip_arr->m, \
                        global_bf[c][b][p].bfchip_arr->k);

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
                bf_free(global_bf[c][b][p].bfchip_arr);

                free(compbuf[c][b][p].comp_arr);
                free(decompbuf[c][b][p].comp_arr);
            }
            free(global_bf[c][b]);

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
            lba = generate_lba('W', pttr);
        }
        else {
            break;
        }
    }
    
    key = lba + offset;
    hashkey = hashing_key(key);
    bf_set(global_bf[chip][blk][offset].bfchip_arr, hashkey);
/*
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

    for(int idx=offset-1; idx>=0; idx--) {
        key = lba + idx;
        hashkey = hashing_key(key);
/*
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
        // Bloomfilter says true
        //if(bf_check((BF*)decompbuf[chip][blk][idx].comp_arr, hashkey) == true) {
        if(bf_check(global_bf[chip][blk][idx].bfchip_arr, hashkey) == true) {
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
    printf("\n### BENCHMARK RESULTS ###\n");
    printf("Total read requests: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("Total false num: %d\n", false_cnt);
    printf("RAF: %.2f\n", (double)(found_cnt + notfound_cnt) / read_cnt);

    uint64_t total_bits = 0;
    printf("Sum of BF assigned bits per chip: ");
    for(int c=0; c<CHIP; c++) {
        sum += bf_man->bits_per_chip[c];
    }
    printf("(SUM: %lu bits, %lu bytes)\n", sum, sum/8);
    total_bits = sum;
    
    sum = 0;
    printf("Sum of BF assigned bits per block in one chip: ");
    for(int b=0; b<BLOCK_PER_CHIP; b++) {
        sum += bf_man->bits_per_blk[0][b];
    }
    printf("(SUM: %lu bits, %lu bytes)\n", sum, sum/8);
    
    sum = 0;
    printf("Sum of BF assigned bits per page in one block: ");
    for(int p=0; p<PAGE_PER_BLOCK; p++) {
        sum += bf_man->bits_per_pg[0][0][p];
    }
    printf("(SUM: %lu bits, %lu bytes)\n", sum, sum/8);

    printf("*** Sum of BF's total assigned: %lu bytes (%lu bits)\n", total_bits/8, total_bits);

    printf("DONE !!\n");
    fflush(stdout);
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
        lba = generate_lba('W', "SEQ");
        bloom_write(lba, val, "SEQ");
        //lba = generate_lba('W', "RAND");
        //bloom_write(lba, val, "RAND");
        
        val++; write_cnt++;
    }

    val = 0;
    while(try_rd--) {
        //lba = generate_lba('R', "SEQ");
        lba = generate_lba('R', "RAND");
        bloom_read(lba);

        val++; read_cnt++;
    }


    FILE* fptr = fopen("origin.bin", "w+"); uint64_t sum=0;
    for(int c=0; c<CHIP; c++) {
        for(int b=0; b<BLOCK_PER_CHIP; b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                fwrite(global_bf[c][b][p].bfchip_arr->body, bf_man->bytes_arr[c][b][p], 1, fptr);
                sum += bf_man->bytes_arr[c][b][p];
            }
        }
    }
    fclose(fptr);
    
    print_stats();
    //printf("%lu\n", sum);

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
