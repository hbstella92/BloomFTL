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

BF** bf_init(int entry, int pg_per_blk);
int bf_set(BF** input, int idx, KEYT key);
bool symbol_check(BF** input, int idx, KEYT key, char* symbol, int symb_length, int front_byte, int front_bit, int start, uint8_t* symb_arr, int symb_arr_sz);
bool bf_check(BF** input, int idx, KEYT key);
void bf_free(BF** input, int pg_per_blk);

uint64_t bf_bits(BF* input);
uint64_t bf_bytes(BF* input);
uint32_t bf_func(BF* input);

BF* bf_cpy(BF*);
void bf_save(BF*);
BF* bf_load();

#endif
