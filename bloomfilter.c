#include "ftl_setting.h"
#include "ftl_data.h"
#include "bloomfilter.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

#define PR_SUCCESS 0.9
#define NUM_CHUNK 10

extern int save_fd;

static inline void BITSET(char *input, char offset){
    (*input) |= (1 << offset);
}

static inline bool BITGET(char input, char offset){
    return input & (1 << offset);
}

// Compressible bf_init (# of hash func is set to 1)
BF** bf_init(int entry, int pg_per_blk) {
    double true_p=0.0, false_p=0.0;
    uint64_t sum_bits=0;

    BF** res = (BF**)malloc(sizeof(BF*) * pg_per_blk);

    for(int p=0; p<pg_per_blk; p++) {
        res[p] = (BF*)malloc(sizeof(BF));
       
        if(p == 0) {
            res[p]->start = 0;
            res[p]->m = 0;
            continue;
        }

        res[p]->n = entry;
        res[p]->k = 1;

        true_p = pow(PR_SUCCESS, (double)1/p);
        false_p = 1 - true_p;
        
        res[p]->p = false_p;
        res[p]->m = ceil(-1 / (log(1-pow(res[p]->p,1/1)) / log(exp(1.0))));
        res[p]->targetsize = res[p]->m / 8;
        if(res[p]->m % 8) {
            res[p]->targetsize++;
        }
        sum_bits += res[p]->m;
        res[p]->start = sum_bits - res[p]->m;
    }
    
    uint64_t sum_bytes = sum_bits / 8;
    if(sum_bits % 8) {
        sum_bytes++;
    }
    res[0]->body = (char*)malloc(sum_bytes);
    memset(res[0]->body, 0, sum_bytes);
    return res;
}

/* Compressible bf_init (variable hash func)
BF* bf_init(int entry, float fpr) {
    if(fpr > 1) { return NULL; }
    BF* res = (BF*)malloc(sizeof(BF));
    res->n = entry;
    res->k = 1;
    res->p = fpr;
    res->m = ceil(-1 / (log(1-pow(res->p,1/1)) / log(exp(1.0))));
    
    int targetsize = res->m / 8;
    //unsigned long long targetsize = res->m / 8;
    if(res->m % 8) {
        targetsize++;
    }
    res->body = (char*)malloc(targetsize);
    memset(res->body, 0, targetsize);
    res->p = fpr;
    res->targetsize = targetsize;
    return res;
}
*/

/* Original bf_init (BF is assigned per page)
BF* bf_init(int entry, float fpr){
	if(fpr>1)
		return NULL;
	BF *res=(BF*)malloc(sizeof(BF));
	res->n=entry;
	res->m=ceil((res->n * log(fpr)) / log(1.0 / (pow(2.0, log(2.0)))));
	res->k=round(log(2.0) * (float)res->m / res->n);
	int targetsize=res->m/8;
	if(res->m%8)
		targetsize++;
	res->body=(char*)malloc(targetsize);
	memset(res->body,0,targetsize);
	res->p=fpr;
	res->targetsize=targetsize;
	return res;
}
*/

/* Modified bf_init (BFs are assigned contiguously)
BF** bf_init(int entry, int pg_per_blk) {
    double true_p=0.0, false_p=0.0;
    uint64_t sum_bits=0;

    BF** res = (BF**)malloc(sizeof(BF*) * pg_per_blk);
    for(int p=0; p<pg_per_blk; p++) {
        res[p] = (BF*)malloc(sizeof(BF));
        if(p == 0) {
            continue;
        }

        res[p]->n = entry;

        true_p = pow(PR_SUCCESS, (double)1/p);
        false_p = 1 - true_p;
        
        res[p]->m = ceil((res[p]->n * log(false_p)) / log(1.0 / (pow(2.0, log(2.0)))));
        res[p]->k = round(log(2.0) * (float)res[p]->m / res[p]->n);
        res[p]->targetsize = res[p]->m / 8;
        if(res[p]->m % 8) {
            res[p]->targetsize++;
        }
        res[p]->p = false_p;
        sum_bits += res[p]->m;
        res[p]->start = sum_bits - res[p]->m;
    }

    uint64_t sum_bytes = sum_bits / 8;
    if(sum_bits % 8) {
        sum_bytes++;
    }
    res[0]->body = (char*)malloc(sum_bytes);
    memset(res[0]->body, 0, sum_bytes);
    return res;
}
*/

/* Modified bf_init (BFs are assigned contiguously)
BF** bf_init(int entry, int pg_per_blk) {
    double true_p=0.0, false_p=0.0;
    int sum_targetbits=0;

    BF** res = (BF**)malloc(sizeof(BF*) * pg_per_blk);
    for(int p=0; p<pg_per_blk; p++) {
        res[p] = (BF*)malloc(sizeof(BF));

        res[p]->n = entry;
        
        if(p == 0) {
            false_p = 0.0001;
        } else {
            true_p = pow(PR_SUCCESS, (double)1/p);
            false_p = 1 - true_p;
        }

        res[p]->m = ceil((res[p]->n * log(false_p)) / log(1.0 / (pow(2.0, log(2.0)))));
        res[p]->k = round(log(2.0) * (float)res[p]->m / res[p]->n);
        int targetsize = res[p]->m / 8;
        if(res[p]->m % 8) {
            targetsize++;
        }

        res[p]->p = false_p;
        res[p]->targetsize = targetsize;
        sum_targetbits += targetsize;
        res[p]->delim = sum_targetbits - res[p]->targetsize;
    }

    res[0]->body = (char*)malloc(sum_targetbits);
    memset(res[0]->body, 0, sum_targetbits);
    return res;
}
*/

void magic_number_set() {
    //
}

void symbol_set(uint64_t bf_bits, int idx, KEYT key, uint8_t* symbol, uint64_t* sym_start, uint64_t* sym_length) {
	int end_byte = (sym_start[idx] + sym_length[idx] - 1) / 8;
	int end_bit = (sym_start[idx] + sym_length[idx] - 1) % 8;
	int symb_arr_sz = end_byte - (sym_start[idx] / 8) + 1;
	int remain_chunk = sym_length[idx];
	uint8_t chunk_sz = sym_length[idx] > end_bit + 1 ? end_bit + 1 : sym_length[idx];
	KEYT h = hashfunction(key) % bf_bits;

    // 1
	if(end_bit == 7) {
		symbol[end_byte] |= h << (8 - chunk_sz);
	}
	else {
		symbol[end_byte] |= h & ((1 << chunk_sz) - 1);
	}

	if(symb_arr_sz == 1) {
		goto task_end;
	}

	h >>= chunk_sz;
	remain_chunk -= chunk_sz;
	chunk_sz = remain_chunk > 8 ? 8 : remain_chunk;

    // 2
	symbol[end_byte - 1] |= h << (8 - chunk_sz);
	if(symb_arr_sz == 2) {
		goto task_end;
	}

	h >>= chunk_sz;
	remain_chunk -= chunk_sz;
	chunk_sz = remain_chunk > 8 ? 8 : remain_chunk;

    // 3
	symbol[end_byte - 2] |= h << (8 - chunk_sz);

task_end:
	return;
}

bool bf_check(BF** input, int idx, KEYT key) {
	if(input[idx] == NULL) return false;
	
    KEYT h;
	int block, offset;
    uint64_t start = input[idx]->start;

	for(uint32_t i=0; i<input[idx]->k; i++){
		//MurmurHash3_x86_32(&key,sizeof(key),i,&h);
		h = hashfunction((key << 19) | (i << 7));
		h %= input[idx]->m;

		block = (start + h) / 8;
		offset = (start + h) % 8;
        
        if(!BITGET(input[0]->body[block], offset)){
			return false;
		}
	}
	return true;
}

void bf_free(BF** input, int pg_per_blk) {
	free(input[0]->body);

    for(int p=0; p<pg_per_blk; p++) {
        free(input[p]);
    }
	free(input);
}

uint64_t bf_bits(BF* input) {
    return input->m;
}

uint64_t bf_bytes(BF* input) {
    uint64_t bytes = input->m / 8;
    if(input->m % 8) {
        bytes++;
    }

	return bytes;
}

uint32_t bf_func(BF* input) {
    return input->k;
}

BF* bf_cpy(BF *src) {
	if(src == NULL) return NULL;

	BF* res=(BF*)malloc(sizeof(BF));
	memcpy(res,src,sizeof(BF));
	res->body=(char *)malloc(res->targetsize);
	memcpy(res->body,src->body,res->targetsize);
	return res;
}

/*
void bf_save(BF* input){
	write(save_fd,input,sizeof(BF));
	write(save_fd,input->body,input->targetsize);
}

BF* bf_load(){
	BF *res=(BF*)malloc(sizeof(BF));
	read(save_fd,res,sizeof(BF));
	res->body=(char*)malloc(res->targetsize);
	read(save_fd,res->body,res->targetsize);
	return res;
}*/
