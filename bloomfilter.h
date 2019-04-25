#ifndef __BLOOM_H__
#define __BLOOM_H__

#include"settings.h"
#include"lsm_settings.h"
#include<stdlib.h>
#include<stdint.h>

typedef struct{
	uint32_t k;
	uint64_t m;
    uint64_t targetsize;
	int n;
	float p;
	char *body;
    uint64_t start; // Start index of BF
}BF;

static inline KEYT hashfunction(KEYT key){
	key ^= key >> 15;
	key *= UINT32_C(0x2c1b3c6d);
	key ^= key >> 12;
	key *= UINT32_C(0x297a2d39);
	key ^= key >> 15;
	/*
	key = ~key + (key << 15); // key = (key << 15) - key - 1;
	key = key ^ (key >> 12);
	key = key + (key << 2);
	key = key ^ (key >> 4);
	key = key * 2057; // key = (key + (key << 3)) + (key << 11);
	key = key ^ (key >> 16);*/
	return key;
}

BF** bf_init(int entry, int pg_per_blk);
int bf_set(BF** input, int idx, KEYT key);
void symbol_set(uint64_t bf_bits, int idx, KEYT key, uint8_t* symbol, uint64_t* sym_start, uint64_t* sym_length, bool seq_flag);

static inline bool symbol_check(uint64_t bf_bits, int idx, KEYT key, uint8_t* symbol, uint64_t* sym_start, uint64_t* sym_length) {
	int end_byte = (sym_start[idx] + sym_length[idx] - 1) / 8;
	int end_bit = (sym_start[idx] + sym_length[idx] - 1) % 8;
	int symb_arr_sz = end_byte - (sym_start[idx] / 8) + 1;
	int remain_chunk = sym_length[idx];
	int chunk_cnt = 0;
	uint8_t chunk_sz = end_bit + 1;
	KEYT h;
	
	h = hashfunction((key << 19)) % bf_bits;

	if(remain_chunk - chunk_sz <= 0){
		chunk_sz = sym_length[idx];
	}

	// 1
	remain_chunk -= chunk_sz;
	if(end_bit == 7){
		if(((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte] >> (8 - chunk_sz)))){
            goto not_exist;
		}
	}
	else{
		if((h ^ symbol[end_byte]) & ((1 << chunk_sz) - 1)){
            goto not_exist;
		}
	}
	chunk_cnt++;
	if(chunk_cnt == symb_arr_sz){
        goto exist;
	}
	h >>= chunk_sz;
	chunk_sz = remain_chunk > 8 ? 8 : remain_chunk;

	// 2
	remain_chunk -= chunk_sz;
	if((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte - 1] >> (8 - chunk_sz))){
        goto not_exist;
	}
	chunk_cnt++;
	if(chunk_cnt == symb_arr_sz){
        goto exist;
	}
	h >>= chunk_sz;
	chunk_sz = remain_chunk > 8 ? 8 : remain_chunk;

	// 3
	if((h & ((1 << chunk_sz) - 1)) ^ (symbol[end_byte - 2] >> (8 - chunk_sz))){
        goto not_exist;
	}

exist:
    return true;

not_exist:
    return false;
}

bool bf_check(BF** input, int idx, KEYT key);
void bf_free(BF** input, int pg_per_blk);

uint64_t bf_bits(BF* input);
uint64_t bf_bytes(BF* input);
uint32_t bf_func(BF* input);

BF* bf_cpy(BF*);
void bf_save(BF*);
BF* bf_load();

#endif
