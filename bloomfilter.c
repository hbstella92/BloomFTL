#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

//#include "ftl_setting.h"
//#include "ftl_data.h"
//#include "ftl_func.h"
#include "bloomfilter.h"

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

extern SBlkManager* sblk_man;
extern GlobalBF** global_bf;
extern BFManager* bf_man;
extern Storage storage;
extern GlobalSymb global_symb;
extern SManager* st_man;
//extern int save_fd;

static inline void BITSET(char *input, char offset){
    (*input) |= (1 << offset);
}

static inline bool BITGET(char input, char offset){
    return input & (1 << offset);
}

void symbol_init() {
	// Allocate SManager
	st_man = (SManager*)malloc(sizeof(SManager));
	st_man->sym_bits_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SBLK);
	st_man->new_pidx_start = (uint64_t*)malloc(sizeof(uint64_t) * TOTAL_SBLK);

	memset(st_man->sym_bits_pg, 0, sizeof(uint64_t) * PAGE_PER_SBLK);

	for(int i=0; i<TOTAL_SBLK; i++) {
		st_man->new_pidx_start[i] = 0;
	}

	st_man->sym_bits_total = st_man->dead_bytes_total = 0;
	st_man->sym_bits_chip = st_man->sym_bits_blk = st_man->sym_bits_super_blk = 0;

	// Allocate symbol delimiter
	st_man->sym_start = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SBLK);
	memset(st_man->sym_start, 0, sizeof(uint64_t) * PAGE_PER_SBLK);

	// Calculate symbol bits
	uint64_t symbol_bit=0, sum=0;

	for(int p=1; p<PAGE_PER_SBLK; p++) {
		symbol_bit = bf_man->bits_per_pg[p];
		symbol_bit = ceil(log(symbol_bit) / log(2));

		if(p == 1){
			symbol_bit += 4;
		}
		st_man->sym_bits_pg[p] = symbol_bit;

#if BIT_CHECK
        printf("page %d\tSYMBOL_bits %lu\n", p, st_man->sym_bits_pg[p]);
#endif

		st_man->sym_start[p] = sum;
		sum += symbol_bit;
	}
#if BIT_CHECK
    printf("\n");
#endif
	
    st_man->sym_bits_super_blk = sum;
	st_man->sym_bits_chip = sum * (SBLK_PER_CHIP);
	st_man->sym_bits_total = st_man->sym_bits_chip * CHIP;

	st_man->sym_bytes_total = sum / 8;
	if(sum % 8) {
		st_man->sym_bytes_total++;
	}

	st_man->dead_bytes_total = PAGE_PER_SBLK / 8;

	// Allocate symbol table
	global_symb.symbol = (uint8_t***)malloc(sizeof(uint8_t**) * CHIP);
	memset(global_symb.symbol, 0, sizeof(uint8_t**) * CHIP);

	for(int c=0; c<CHIP; c++) {
		global_symb.symbol[c] = (uint8_t**)malloc(sizeof(uint8_t*) * SBLK_PER_CHIP);
		memset(global_symb.symbol[c], 0, sizeof(char*) * SBLK_PER_CHIP);

		for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
			global_symb.symbol[c][sb] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
			memset(global_symb.symbol[c][sb], 0, sizeof(uint8_t) * st_man->sym_bytes_total);
		}
	}

	global_symb.lba_flag = (int**)malloc(sizeof(int*) * TOTAL_SBLK);
	global_symb.ppa_flag = (int**)malloc(sizeof(int*) * TOTAL_SBLK);
	for(int sb=0; sb<TOTAL_SBLK; sb++) {
		global_symb.lba_flag[sb] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
		global_symb.ppa_flag[sb] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);

		memset(global_symb.lba_flag[sb], 0, sizeof(int) * PAGE_PER_SBLK);
		memset(global_symb.ppa_flag[sb], 0, sizeof(int) * PAGE_PER_SBLK);
	}
}

void symbol_destroy() {
	free(st_man->sym_bits_pg);
	free(st_man->new_pidx_start);
	free(st_man->sym_start);
	free(st_man);

	for(int c=0; c<CHIP; c++) {
		for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
			free(global_symb.symbol[c][sb]);
		}

		free(global_symb.symbol[c]);
	}

	for(int sb=0; sb<TOTAL_SBLK; sb++) {
		free(global_symb.lba_flag[sb]);
		free(global_symb.ppa_flag[sb]);
	}

	free(global_symb.symbol);
	free(global_symb.lba_flag);
	free(global_symb.ppa_flag);
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

void symbol_resymbolize(uint32_t pbn, int chip, int way, int chnl, int blk, int valid_start, int num_flush, int sb) {
	uint32_t hashkey;
	int superblk = chip * SBLK_PER_CHIP + sb;

	storage.chip_arr[way][chnl].empty[sb] += num_flush;

	free(global_symb.symbol[chip][sb]);
	sblk_man->num_bf[superblk] = 0;

	global_symb.symbol[chip][sb] = (uint8_t*)malloc(sizeof(uint8_t) * st_man->sym_bytes_total);
	memset(global_symb.symbol[chip][sb], 0, sizeof(uint8_t) * st_man->sym_bytes_total);

	int ppa_end_idx = storage.chip_arr[way][chnl].empty[sb];
	int new_bf_idx = 1;
    int rb_chunk = 0;

	// Set new symbol and Set dead bit for LBA and PPA
	for(int pa=0, b=0; pa<ppa_end_idx, b<SBLK_SZ; pa++) {
		bool make_new_bf = false;

		int local_pa = pa % PAGE_PER_BLOCK;

		if(pa < valid_start) {
			if(global_symb.ppa_flag[superblk][pa] == 1) {
				make_new_bf = true;
			}
		}
		else if(pa > valid_start) {
			if(local_pa == 0) {
				if(storage.chip_arr[way][chnl].data_blk[sb][b-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
						storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[0].oob) {
					global_symb.ppa_flag[superblk][pa] = 0;
                    rb_chunk++;
				}
				else {
					make_new_bf = true;
				}
			}
			else {
				if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 ==
						storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
					global_symb.ppa_flag[superblk][pa] = 0;
                    rb_chunk++;
				}
				else {
					make_new_bf = true;
				}
			}
		}
		else {
			if(pa != 0) {
				if(global_symb.ppa_flag[superblk][pa-1] == 1) {
					if(local_pa == 0) {
						if(storage.chip_arr[way][chnl].data_blk[sb][b-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
								storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[0].oob) {
							global_symb.ppa_flag[superblk][pa] = 0;
                            rb_chunk++;
						}
						else {
							make_new_bf = true;
						}
					}
					else {
						if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 == 
								storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
							global_symb.ppa_flag[superblk][pa] = 0;
                            rb_chunk++;
						}
						else {
							make_new_bf = true;
						}
					}
				}
				else {
					make_new_bf = true;
				}
			}
		}

        if(rb_chunk == MAX_RB) {
            make_new_bf = true;
        }

		if(make_new_bf == true && pa != 0) {
			hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob);
			symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
					global_symb.symbol[chip][sb], st_man->sym_start, st_man->sym_bits_pg);

			global_symb.ppa_flag[superblk][pa] = 1;

			new_bf_idx++;

            rb_chunk = 0;
		}

		if(!((pa + 1) % PAGE_PER_BLOCK)) {
			b++;
		}

		if(pa + 1 == storage.chip_arr[way][chnl].empty[sb]) {
			break;
		}
	}

	sblk_man->num_bf[superblk] = new_bf_idx - 1;
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

#if SYMMETRIC
        res[p]->p = (double)0.1*2/(PAGE_PER_SBLK + 1);
#else
        true_p = pow(PR_SUCCESS, (double)1/p);
        false_p = 1 - true_p;
        
        res[p]->p = false_p;
#endif

        res[p]->m = ceil(-1 / (log(1-pow(res[p]->p,1/1)) / log(exp(1.0))));

#if BIT_CHECK
        printf("page %d\tBF_Bits %lu\n", p, res[p]->m);
#endif

        res[p]->targetsize = res[p]->m / 8;
        if(res[p]->m % 8) {
            res[p]->targetsize++;
        }
        sum_bits += res[p]->m;
        res[p]->start = sum_bits - res[p]->m;
    }
#if BIT_CHECK
    printf("\n");
#endif

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

/*
BF* bf_cpy(BF *src) {
	if(src == NULL) return NULL;

	BF* res=(BF*)malloc(sizeof(BF));
	memcpy(res,src,sizeof(BF));
	res->body=(char *)malloc(res->targetsize);
	memcpy(res->body,src->body,res->targetsize);
	return res;
}

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
