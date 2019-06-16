#ifndef __BLOOM__
#define __BLOOM__

#include <stdlib.h>
#include <stdint.h>

#include "ftl_setting.h"
#include "ftl_data.h"
#include "ftl_func.h"

static inline KEYT hashfunction(KEYT key) {
    key ^= key >> 15;
	key *= 2246822519U;
	key ^= key >> 13;
	key *= 3266489917U;
	key ^= key >> 16;
	
    return key;
}

void symbol_init();
void symbol_destroy();
void symbol_set(uint64_t bf_bits, int idx, KEYT key, uint8_t* symbol, \
        uint64_t* sym_start, uint64_t* sym_length);

static inline bool symbol_check(uint64_t bf_bits, KEYT key, uint8_t* symbol, \
        uint64_t sym_start, uint64_t sym_length) {
    int end_byte = (sym_start + sym_length - 1) / 8;
    int end_bit = (sym_start + sym_length - 1) % 8;
    int symb_arr_sz = end_byte - (sym_start / 8) + 1;
    uint8_t chunk_sz = sym_length > end_bit + 1 ? end_bit + 1 : sym_length;
    KEYT h = hashfunction(key) % bf_bits;

    // 1
	if(end_bit == 7) {
        if(((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte] >> (8 - chunk_sz)))) {
			goto not_exist;
		}
	}
    else{
        if((h ^ symbol[end_byte]) & ((1 << chunk_sz) - 1)){
            goto not_exist;
        }
    }

	if(symb_arr_sz == 1) {
		goto exist;
	}

	end_byte--;
	h >>= chunk_sz;
	sym_length -= chunk_sz;
	chunk_sz = sym_length > 8 ? 8 : sym_length;

    // 2
    if((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte] >> (8 - chunk_sz))){
        goto not_exist;
    }

	if(symb_arr_sz == 2) {
		goto exist;
	}

	end_byte--;
	h >>= chunk_sz;
	sym_length -= chunk_sz;
	chunk_sz = sym_length > 8 ? 8 : sym_length;

    // 3
    if((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte] >> (8 - chunk_sz))){
        goto not_exist;
    }

exist:
    return true;

not_exist:
    return false;
}

void symbol_resymbolize(uint32_t pbn, int chip, int way, int chnl, int blk, \
        int valid_start, int num_flush, int sb);

BF** bf_init(int entry, int pg_per_blk);
bool bf_check(BF** input, int idx, KEYT key);
void bf_free(BF** input, int pg_per_blk);

uint64_t bf_bits(BF* input);
uint64_t bf_bytes(BF* input);
uint32_t bf_func(BF* input);

// Fibonacci hash
static inline uint32_t hashing_key(uint32_t key) {
    return (uint32_t)((0.618033887 * key) * 1024);
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

/*
BF* bf_cpy(BF*);
void bf_save(BF*);
BF* bf_load();
*/

#endif
