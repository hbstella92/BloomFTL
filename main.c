#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "settings.h"
#include "lsm_settings.h"
#include "bloomfilter.h"
#include "sha1.h"

#define CHANNEL 8
#define WAY 8
#define NUM_PAGES 128
#define NUM_BLOCKS (CHANNEL*WAY)
#define RANGE_LBA (NUM_BLOCKS*NUM_PAGES)

uint32_t data;
int read_cnt;
int write_cnt;
int read_cache_hit;
int read_cache_miss;

int true_cnt;
int false_cnt;

int found_cnt;
int notfound_cnt;

BF** global_bf;

typedef struct block {
    uint32_t* page;
    uint32_t* oob;
    uint32_t empty;
} Block;

Block data_block[CHANNEL][WAY];

uint32_t generate_lba(char* mode) {
    if(!strcmp(mode, "SEQ")) {
        return data;
    }
    else if(!strcmp(mode, "RAND")) {
        return rand() % RANGE_LBA;
    }
}

uint32_t* hashing_key(uint32_t key) {
    uint32_t* hashkey;
    char* string = (char*)&key;
    char result[21];
    char hexresult[41];
    size_t offset;

    SHA1(result, string, strlen(string));

    for(offset=0; offset<20; offset++) {
        sprintf((hexresult + (2*offset)), "%02x", result[offset] & 0xff);
    }

    hashkey = (uint32_t*)hexresult;
    return hashkey;
}

void bloom_write(uint32_t lba, uint32_t data) {
    uint32_t pbn, key, *hashkey;
    int way, chnl, empty;

    pbn = (lba & ((1 << 6) - 1));
    way = pbn / CHANNEL;
    chnl = pbn % CHANNEL;
    empty = data_block[way][chnl].empty;
    // SEQWR
    key = lba + empty;
    hashkey = hashing_key(key);

    data_block[way][chnl].page[empty] = data;
    data_block[way][chnl].oob[empty] = lba;
    data_block[way][chnl].empty++;
  
    bf_set(global_bf[pbn], *hashkey);
}

int bloom_read(uint32_t lba) {
    uint32_t pbn, key, data, *hashkey;
    int way, chnl, empty;

    pbn = (lba & ((1 << 6) - 1));
    way = pbn / CHANNEL;
    chnl = pbn % CHANNEL;
    empty = data_block[way][chnl].empty;

    for(int idx=empty-1; idx>=0; idx--) {
        key = lba + idx;
        hashkey = hashing_key(key);
        
        if(bf_check(global_bf[pbn], *hashkey) == true) { // Bloomfilter true - true or false
            if(data_block[way][chnl].oob[idx] == lba) { // really true
                found_cnt++;
                return data;
            }
            else { // really false
                notfound_cnt++;
                continue;
            }
        }
    }
}

int main() {
    int try_read, try_write;
    uint32_t lba, pbn;

    try_read = try_write = RANGE_LBA;
    read_cnt = write_cnt = read_cache_hit = read_cache_miss = 0;
    true_cnt = false_cnt = 0;
    found_cnt = notfound_cnt = 0;

    global_bf = (BF**)malloc(sizeof(BF*) * NUM_BLOCKS);
    for(int i=0; i<NUM_BLOCKS; i++) {
        global_bf[i] = bf_init(NUM_PAGES, 0.1);
    }

    for(int i=0; i<WAY; i++) {
        for(int j=0; j<CHANNEL; j++) {
            data_block[i][j].page = (uint32_t*)malloc(sizeof(uint32_t) * NUM_PAGES);
            data_block[i][j].oob = (uint32_t*)malloc(sizeof(uint32_t) * NUM_PAGES);
        }
    }

    for(int i=0; i<WAY; i++) {
        for(int j=0; j<CHANNEL; j++) {
            for(int k=0; k<NUM_PAGES; k++) {
                data_block[i][j].page[k] = data_block[i][j].oob[k] = -1;
            }
            data_block[i][j].empty = 0;
        }
    }

    data = 0;
    while(try_write--) {
        lba = generate_lba("SEQ");
        
        bloom_write(lba, data);
        data++;
        write_cnt++;
    }
/*
    printf("TEST write\n");
    for(int i=0; i<WAY; i++) {
        for(int j=0; j<CHANNEL; j++) {
            printf("In block %d\n", i * CHANNEL + j);
            for(int k=0; k<NUM_PAGES; k++) {
                printf("%d\t", data_block[i][j].page[k]);
            }
            printf("\n");
        }
    }
    printf("DONE !!\n\n");
*/

    srand((unsigned int)time(NULL));

    while(try_read--) {
        lba = generate_lba("RAND");
        
        bloom_read(lba);
        read_cnt++;
    }

    printf("\nTEST\n");
    printf("NUM READ: %d\n", read_cnt);
    printf("Total found num: %d\n", found_cnt);
    printf("Total not-found num: %d\n", notfound_cnt);
    printf("RAF: %.2f\n", (float)(found_cnt + notfound_cnt) / read_cnt);
    printf("DONE !!\n");

    for(int i=0; i<WAY; i++) {
        for(int j=0; j<CHANNEL; j++) {
            free(data_block[i][j].page);
            free(data_block[i][j].oob);
        }
    }

    for(int i=0; i<NUM_BLOCKS; i++) {
        bf_free(global_bf[i]);
    }
    free(global_bf);

    return 0;
}
