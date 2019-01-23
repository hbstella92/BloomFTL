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
#define NUM_PAGES 128
#define NUM_BLOCKS (CHANNEL*WAY)
#define RANGE_LBA (NUM_BLOCKS*NUM_PAGES)

typedef struct check {
    uint32_t lba;
    uint32_t d;
} check;

check written[RANGE_LBA];
check reading[RANGE_LBA];

uint32_t lba_arr[RANGE_LBA];
uint32_t data;

int num;
int read_cnt;
int write_cnt;

int true_cnt;
int false_cnt;
int found_cnt;
int notfound_cnt;

BF** global_bf;

typedef struct block {
    uint32_t page[NUM_PAGES];
    uint32_t oob[NUM_PAGES];
    int empty;
} Block;

Block** data_block;

uint32_t generate_lba(char* req, char* mode) {
    if(!strcmp(req, "SEQ")) {
        if(!strcmp(mode, "WR")) {
            return data;
        }
        return lba_arr[data];
    }
    else if(!strcmp(req, "RAND")) {
        if(!strcmp(mode, "RD")) {
            int rand_idx = rand() % RANGE_LBA;
            return lba_arr[rand_idx];
        }
        return rand() % RANGE_LBA;
    }
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

void print_stats() {
    /*
    printf("\nWRITE REQUEST\n");
    for(int i=0; i<RANGE_LBA; i++) {
        printf("%d\t%u\n", i, wr_hash_arr[i]);
    }

    printf("\nREAD REQUEST\n");
    for(int i=0; i<RANGE_LBA; i++) {
        printf("%d\t%u\n", i, rd_hash_arr[i]);
    }
    */
}

void bloom_write(char* req, uint32_t data) {
    uint32_t lba, pbn, key, hashkey;
    int way, chnl, empty;

    do {
        lba = generate_lba(req, "WR");
        pbn = (lba & ((1 << 6) - 1));
        way = pbn / CHANNEL;
        chnl = pbn % CHANNEL;
        empty = data_block[way][chnl].empty;
    } while(empty >= NUM_PAGES);

    // SEQWR
    //key = hashing_key(lba) + empty;
    key = lba + empty;
    hashkey = hashing_key(key);
    //wr_hash_arr[write_cnt] = hashkey;

    //bf_set(global_bf[pbn], hashkey);
    bf_set(global_bf[pbn*NUM_PAGES+empty], hashkey);
    data_block[way][chnl].page[empty] = data;
    data_block[way][chnl].oob[empty] = lba;
written[num].lba = lba;
written[num].d = data;
    data_block[way][chnl].empty++;
    lba_arr[num++] = lba;
}

void bloom_read(char* req) {
    uint32_t lba, pbn, key, data, hashkey;
    int way, chnl, empty;

    lba = generate_lba(req, "RD");
    pbn = (lba & ((1 << 6) - 1));
    way = pbn / CHANNEL;
    chnl = pbn % CHANNEL;
    empty = data_block[way][chnl].empty;

    for(int idx=empty-1; idx>=0; idx--) {
        //key = hashing_key(lba) + idx;
        key = lba + idx;
        hashkey = hashing_key(key);
        //rd_hash_arr[read_cnt] = hashkey;
        
        if(bf_check(global_bf[pbn*NUM_PAGES+idx], hashkey) == true) { // Bloomfilter true - true or false
        //if(bf_check(global_bf[pbn], hashkey) == true) { // Bloomfilter true - true or false
            if(data_block[way][chnl].oob[idx] == lba) { // really true
reading[found_cnt].lba = lba;
reading[found_cnt].d = data_block[way][chnl].page[idx];
                found_cnt++;
                return;
            }
            else { // really false
                notfound_cnt++;
                continue;
            }
        }
        else {
            false_cnt++;
        }
    }
}

int main() {
    int try_read, try_write;
    uint32_t lba, pbn;
    double ptv_p=0.0, false_p=0.0;
    uint64_t bytes=0, sum_bytes=0;

    try_read = try_write = RANGE_LBA;
    num = read_cnt = write_cnt = 0;
    true_cnt = false_cnt = 0;
    found_cnt = notfound_cnt = 0;

    global_bf = (BF**)malloc(sizeof(BF*) * RANGE_LBA);
    for(int b=0; b<NUM_BLOCKS; b++) {
        for(int p=0; p<NUM_PAGES; p++) {
            if((p == 0) && (b == 0)) {
                global_bf[0] = bf_init(1, 0.001);
            }
            else {
                ptv_p = pow(0.9, (double)1/p);
                false_p = 1 - ptv_p;
                
                global_bf[b*NUM_PAGES+p] = bf_init(1, false_p);
                bytes += bf_bits(1, false_p);
            }
        }
        
//        printf("Block %d. bloom bytes: %lu\n", b, bytes);
        sum_bytes += bytes;
        bytes=0;
    }

    /*
     * TRADITIONAL BF
     *
    global_bf = (BF**)malloc(sizeof(BF*) * NUM_BLOCKS);
    for(int i=0; i<NUM_BLOCKS; i++) {
        global_bf[i] = bf_init(NUM_PAGES, 0.1);
    }
    */

    data_block = (Block**)malloc(sizeof(Block*) * WAY);
    for(int i=0; i<WAY; i++) {
        data_block[i] = (Block*)malloc(sizeof(Block) * CHANNEL);
    }
    
    for(int i=0; i<WAY; i++) {
        for(int j=0; j<CHANNEL; j++) {
            for(int k=0; k<NUM_PAGES; k++) {
                data_block[i][j].page[k] = data_block[i][j].oob[k] = 9999;
            }
            data_block[i][j].empty = 0;
        }
    }

    srand((unsigned int)time(NULL));

    data = 0;
    while(try_write--) {
        bloom_write("RAND", data);
        data++;
        write_cnt++;
    }

    data = 0;
    while(try_read--) {
        bloom_read("SEQ");
        data++;
        read_cnt++;
    }

    int success=0;
    printf("\n\n\n");
    for(int i=RANGE_LBA; i>=0; i--) {
        for(int j=RANGE_LBA; j>=0; j--) {
            if((reading[i].lba == written[j].lba) && (reading[i].d == written[j].d)) {
                success++;
            }
        }
    }
    printf("SUCCESS: %d\n", success);
    printf("\n\n\n");

    printf("TEST\n");
    printf("NUM READ: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("Total error num: %d\n", false_cnt);
    printf("RAF: %.2f\n", (float)(found_cnt + notfound_cnt) / read_cnt);
    printf("Sum of bloom bytes: %lu\n", sum_bytes);
    printf("DONE !!\n");

    for(int i=0; i<WAY; i++) {
        free(data_block[i]);
    }
    free(data_block);

    for(int i=0; i<RANGE_LBA; i++) {
        bf_free(global_bf[i]);
    }
    free(global_bf);

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
