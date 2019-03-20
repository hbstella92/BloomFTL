#ifndef __BLOOM_H__
#define __BLOOM_H__
#include"settings.h"
#include"lsm_settings.h"
#include<stdlib.h>
#include<stdint.h>

typedef struct{
	uint32_t k;
	int m;
	int targetsize;
	int n;
	float p;
	char *body;
    int delim;
}BF;

BF** bf_init(int, int);
uint32_t bf_func(BF*);
void bf_free(BF**, int);
uint64_t bf_bits(BF*);
uint64_t bf_bytes(BF*);
void bf_set(BF**, int, KEYT);
bool bf_check(BF**, int, KEYT);

BF* bf_cpy(BF*);
void bf_save(BF*);
BF* bf_load();
#endif
