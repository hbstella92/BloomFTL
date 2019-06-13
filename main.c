#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "ftl_setting.h"
#include "ftl_func.h"
#include "bloomfilter.h"

#if ARM
#include "sha256-arm.h"
#else
#include "sha256.h"
#endif

SBlkManager* sblk_man;
GlobalBF** global_bf;
BFManager* bf_man;
Storage storage;
GlobalSymb global_symb;
SManager* st_man;
#if BUFFERED
BlockBuffer** block_buffer;
#endif

double read_check=0.0;
int read_loop=0;
double w_time=0.0;
double r_time=0.0;

uint32_t val;

uint32_t mask;

int read_cnt;
int write_cnt;
int true_cnt;
long int false_cnt;
int found_cnt;
int notfound_cnt;
int evict_blk_cnt;
int disk_write_cnt;

void bloom_init() {
	double true_p=0.0, false_p=0.0;

	mask = (int)(PAGE_PER_SBLK*(1 - OP));

	// SBlkManager
	sblk_man = (SBlkManager*)malloc(sizeof(SBlkManager));
	sblk_man->num_bf = (int*)malloc(sizeof(int) * TOTAL_SBLK);
	for(int i=0; i<TOTAL_SBLK; i++) {
		sblk_man->num_bf[i] = 0;
	}

	// BFManager
	bf_man = (BFManager*)malloc(sizeof(BFManager));
	bf_man->bits_per_pg = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SBLK);
	bf_man->bytes_arr = (uint64_t*)malloc(sizeof(uint64_t) * PAGE_PER_SBLK);

	for(int p=0; p<PAGE_PER_SBLK; p++) {
		bf_man->bits_per_pg[p] = 0;
		bf_man->bytes_arr[p] = 0;
	}

	// BF
	global_bf = (GlobalBF**)malloc(sizeof(GlobalBF*) * CHIP);
	for(int c=0; c<CHIP; c++) {
		global_bf[c] = (GlobalBF*)malloc(sizeof(GlobalBF) * SBLK_PER_CHIP);

		for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
			global_bf[c][sb].bfchip_arr = bf_init(1, PAGE_PER_SBLK);
		}
	}

	for(int p=0; p<PAGE_PER_SBLK; p++) {
		if(p != 0) {
			bf_man->bits_per_pg[p] = bf_bits(global_bf[0][0].bfchip_arr[p]);

			int targetsize = bf_man->bits_per_pg[p] / 8;
			if(bf_man->bits_per_pg[p] % 8) {
				targetsize++;
			}
			bf_man->bytes_arr[p] = targetsize;
		}
	}

	// All blocks and pages
	storage.chip_arr = (Chip**)malloc(sizeof(Chip*) * WAY);
	for(int w=0; w<WAY; w++) {
		storage.chip_arr[w] = (Chip*)malloc(sizeof(Chip) * CHANNEL);

		for(int ch=0; ch<CHANNEL; ch++) {
			storage.chip_arr[w][ch].data_blk = (Block**)malloc(sizeof(Block*) * SBLK_PER_CHIP);
			storage.chip_arr[w][ch].empty = (int*)malloc(sizeof(int) * SBLK_PER_CHIP);

			for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
				storage.chip_arr[w][ch].data_blk[sb] = (Block*)malloc(sizeof(Block) * SBLK_SZ);

				storage.chip_arr[w][ch].empty[sb] = 0;

				for(int b=0; b<SBLK_SZ; b++) {
					storage.chip_arr[w][ch].data_blk[sb][b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

					memset(storage.chip_arr[w][ch].data_blk[sb][b].page_arr, 0, sizeof(Page) * PAGE_PER_BLOCK);
				}
			}
		}
	}

#if BUFFERED
	block_buffer = (BlockBuffer**)malloc(sizeof(BlockBuffer*) * CHIP);
	for(int c=0; c<CHIP; c++) {
		block_buffer[c] = (BlockBuffer*)malloc(sizeof(BlockBuffer) * BLOCK_PER_CHIP);

		for(int b=0; b<BLOCK_PER_CHIP; b++) {
			block_buffer[c][b].lba = (uint32_t*)malloc(sizeof(uint32_t) * BUF_SZ);
			memset(block_buffer[c][b].lba, 0, sizeof(uint32_t) * BUF_SZ);

			block_buffer[c][b].value = (uint32_t*)malloc(sizeof(uint32_t) * BUF_SZ);
			memset(block_buffer[c][b].value, 0, sizeof(uint32_t) * BUF_SZ);

			block_buffer[c][b].buf_sz = 0;
		}
	}
#endif
}

void bloom_destroy() {
	free(sblk_man->num_bf);
	free(sblk_man);

	for(int c=0; c<CHIP; c++) {
		for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
			bf_free(global_bf[c][sb].bfchip_arr, PAGE_PER_SBLK);
#if BUFFERED
			free(block_buffer[c][b].lba);
			free(block_buffer[c][b].value);
#endif
		}

		free(global_bf[c]);
#if BUFFERED
		free(block_buffer[c]);
#endif
	}

	free(bf_man->bytes_arr);
	free(bf_man->bits_per_pg);
	free(bf_man);
	free(global_bf);
#if BUFFERED
	free(block_buffer);
#endif

	for(int w=0; w<WAY; w++) {
		for(int ch=0; ch<CHANNEL; ch++) {
			for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
				for(int b=0; b<SBLK_SZ; b++) {
					free(storage.chip_arr[w][ch].data_blk[sb][b].page_arr);
				}

				free(storage.chip_arr[w][ch].data_blk[sb]);
			}

			free(storage.chip_arr[w][ch].data_blk);
			free(storage.chip_arr[w][ch].empty);
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

/* 
 *
	NOTE: If SEQ_FLAG is delivered, these flags does not work!
	W_UNIQUE: Random && Unique
	R_UNIQUE: Random && Unique (If W_UNIQUE is not set, R_UNIQUE doesn't work)
	R_RD_ORDER:
		ex) If W_UNIQUE is set and writes LBA 3-1-2, R_RD_ORDER reads LBA 2-1-3
	R_RD_GEN: Random && Duplicated
 */
void make_test_set(uint32_t *w_arr, uint32_t *r_arr, char *w_t, char *r_t, uint8_t options, uint32_t test_size) {
	uint8_t *page_usage;
	uint32_t make_cnt, lba, pbn;
	uint8_t read = 1, write = 1;

	if(!strcmp(w_t, "RAND")) {
		write = 0;
	}

	// Generate write LBA
	if(!write && !(options & W_UNIQUE)) {
		// Check whether the physical block corresponding to LBA is full or not (Because there is no GC algorithm)
		make_cnt = 0;
		page_usage = (uint8_t*)malloc(TOTAL_PAGE);
		memset(page_usage, 0, TOTAL_PAGE);

		while(make_cnt < test_size) {
			lba = rand() % TOTAL_PAGE;
			pbn = lba / mask; // random write

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
	if(!strcmp(r_t, "RAND")) {
		read = 0;
	}

	// Generate read LBA
	if(read) { // Seq read
		memcpy(r_arr, w_arr, sizeof(uint32_t) * test_size);
		if(!write) {
			qsort(r_arr, test_size, sizeof(uint32_t), compare);
		}
	}
	else { // Rand read
		if((write && !(options & R_UNIQUE)) ||
				(!write && ((options & (W_UNIQUE | R_UNIQUE)) == W_UNIQUE ||\
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

// For debugging
void print_byte_as_bit(unsigned char* val, size_t bytes) {
	for(int i=0; i<bytes; i++) {
		for(int j=7; j>=0; j--) {
			printf("%c", (val[i] & (1 << j)) ? '1' : '0');
		} printf(" ");
	} printf("\n");
}

int lba_compare(const void* a, const void* b) {
	Page* one = (Page*)a;
	Page* two = (Page*)b;

	if(one->oob < two->oob) {
		return -1;
	}
	else if(one->oob > two->oob) {
		return 1;
	}
	else {
		return 0;
	}
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
				}
				else {
					make_new_bf = true;
				}
			}
			else {
				if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 ==
						storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
					global_symb.ppa_flag[superblk][pa] = 0;
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
						}
						else {
							make_new_bf = true;
						}
					}
					else {
						if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa-1].oob + 1 == 
								storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob) {
							global_symb.ppa_flag[superblk][pa] = 0;
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

		if(make_new_bf == true && pa != 0) {
			hashkey = hashing_key(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[local_pa].oob);
			symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, hashkey + new_bf_idx, \
					global_symb.symbol[chip][sb], st_man->sym_start, st_man->sym_bits_pg);

			global_symb.ppa_flag[superblk][pa] = 1;

			new_bf_idx++;
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

void bloom_gc_revision(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, Page* rebloom_list, int num_rebloom, int* sb) {
    int superblk = (*chip) * SBLK_PER_CHIP + (*sb);
    Page *candi_list, *valid_list, prev_page;
    int num_valid = 0, candi_sz = 0;
    int gidx = 0;

    if(!rebloom_list) {
        // Find valid pages in evict block
        candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
        valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

        memset(candi_list, 0, sizeof(Page) * PAGE_PER_BLOCK);
        memset(valid_list, 0, sizeof(Page) * PAGE_PER_BLOCK);

        for(int i=0; i<PAGE_PER_BLOCK; i++) {
            if(global_symb.ppa_flag[superblk][i] == 1) {
                candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i];
            }
            else if((i != 0) && (storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i-1].oob + 1 ==
                    storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i].oob)) {
                candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i];
            }
        }

        qsort(candi_list, candi_sz, sizeof(Page), lba_compare);

        prev_page = valid_list[num_valid++] = candi_list[0];
        for(int i=1; i<candi_sz; i++) {
            if(prev_page.oob != candi_list[i].oob) {
                valid_list[num_valid++] = candi_list[i];
            }

            prev_page = candi_list[i];
        }

        // Invalid with ppa_flag
        for(int b=1; b<SBLK_SZ; b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                gidx = b * PAGE_PER_BLOCK + p;

                if(global_symb.ppa_flag[superblk][gidx] == 1) {
                    for(int i=0; i<num_valid; i++) {
                        if(valid_list[i].oob == storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                            global_symb.ppa_flag[superblk][gidx] = 0;
                            break;
                        }
                    }
                }
            }
        }

#if DEBUG
        printf("CANDI LIST\n");
        for(int i=0; i<candi_sz; i++) {
            printf("[%d]\tlba: %u\n", i, candi_list[i].oob);
        }
        printf("\n");
        printf("VALID LIST\n");
        for(int i=0; i<num_valid; i++) {
            printf("[%d]\tlba: %u\n", i, valid_list[i].oob);
        }
        printf("\n");
#endif

        free(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr);
        evict_blk_cnt++;
        
        // Re-allocate new block and Copy valid pages
        storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
        memset(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr, 0, sizeof(Page) * PAGE_PER_BLOCK);
        storage.chip_arr[*way][*chnl].empty[*sb] = (3 * PAGE_PER_BLOCK);

        Page* tmp = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
        memset(tmp, 0, sizeof(Page) * PAGE_PER_BLOCK);
       
        tmp = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;
        for(int b=0; b<SBLK_SZ-1; b++) {
            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
                storage.chip_arr[*way][*chnl].data_blk[*sb][b+1].page_arr;
        }
        storage.chip_arr[*way][*chnl].data_blk[*sb][SBLK_SZ-1].page_arr = tmp;

        memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][SBLK_SZ-1].page_arr, valid_list, sizeof(Page) * num_valid);

        // Manage ppa_flag
        int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
        memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SBLK);

        memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK], sizeof(int) * (PAGE_PER_BLOCK * (SBLK_SZ - 1)));
        free(global_symb.ppa_flag[superblk]);

        global_symb.ppa_flag[superblk] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
        memset(global_symb.ppa_flag[superblk], 0, sizeof(int) * PAGE_PER_SBLK);
        memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SBLK);

        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, PAGE_PER_BLOCK * 3, num_valid, *sb);

        disk_write_cnt += num_valid;

        free(tmp_ppa_flag);
        //free(tmp);
        free(valid_list);
        free(candi_list);
    }
    else {
        candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SBLK);
        valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SBLK);

        memset(candi_list, 0, sizeof(Page) * PAGE_PER_SBLK);

        memcpy(candi_list, rebloom_list, sizeof(Page) * num_rebloom);
        candi_sz += num_rebloom;

        int evict_blk = -1;
        
        do {
            if(evict_blk != -1) {
                memset(candi_list, 0, sizeof(Page) * PAGE_PER_SBLK);
                memcpy(candi_list, valid_list, sizeof(Page) * num_valid);
                candi_sz = num_valid;
            }

            memset(valid_list, 0, sizeof(Page) * PAGE_PER_SBLK);
            num_valid = 0;
            
            evict_blk++;

            for(int i=0; i<PAGE_PER_BLOCK; i++) {
                gidx = evict_blk * PAGE_PER_BLOCK + i;

                if(gidx == (*ppn)) {
                    break;
                }
            
                if(global_symb.ppa_flag[superblk][gidx] == 1) {
                    candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i];
                }
                else if(i != 0) {
                    if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i-1].oob + 1 ==
                            storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i].oob) {
                        candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i];
                    }
                }
                else if(evict_blk != 0) {
                    if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
                            storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i].oob) {
                        candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[i];
                    }
                }
            }

            qsort(candi_list, candi_sz, sizeof(Page), lba_compare);

            prev_page = valid_list[num_valid++] = candi_list[0];
            for(int i=1; i<candi_sz; i++) {
                if(prev_page.oob != candi_list[i].oob) {
                    valid_list[num_valid++] = candi_list[i];
                }

                prev_page = candi_list[i];
            }

            for(int b=evict_blk+1; b<SBLK_SZ; b++) {
                for(int p=0; p<PAGE_PER_BLOCK; p++) {
                    gidx = b * PAGE_PER_BLOCK + p;

                    if(global_symb.ppa_flag[superblk][gidx] == 1) {
                        for(int i=0; i<num_valid; i++) {
                            if(valid_list[i].oob == storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                                global_symb.ppa_flag[superblk][gidx] = 0;
                                break;
                            }
                        }
                    }
                }
            }
        } while(((SBLK_SZ - 1 - evict_blk) * PAGE_PER_BLOCK) + num_valid > PAGE_PER_SBLK);

#if DEBUG
        printf("VALID LIST\n");
        for(int i=0; i<num_valid; i++) {
            printf("[%d]\tlba: %u\n", i, valid_list[i].oob);
        }
        printf("\n");
#endif

        evict_blk++;
        
        for(int b=0; b<evict_blk; b++) {
            free(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr);

            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
            memset(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr, 0, sizeof(Page) * PAGE_PER_BLOCK);
        }
        
        evict_blk_cnt += evict_blk;

        storage.chip_arr[*way][*chnl].empty[*sb] -= (PAGE_PER_BLOCK * evict_blk);

        Page* tmp = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
        memset(tmp, 0, sizeof(Page) * PAGE_PER_BLOCK);

        for(int i=0; i<evict_blk; i++) {
            tmp = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;

            for(int b=0; b<SBLK_SZ-1; b++) {
                storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
                    storage.chip_arr[*way][*chnl].data_blk[*sb][b+1].page_arr;
            }
            storage.chip_arr[*way][*chnl].data_blk[*sb][SBLK_SZ-1].page_arr = tmp;
        }

        // Flush
        int remain_valid = num_valid;
        int start_blk = storage.chip_arr[*way][*chnl].empty[*sb] / PAGE_PER_BLOCK;
        int start_idx = storage.chip_arr[*way][*chnl].empty[*sb] % PAGE_PER_BLOCK;
        int j = 0, v = 0;
        int dat = 0;

        while(remain_valid > 0) {
            if(j == 0) {
                dat = PAGE_PER_BLOCK - start_idx;

                memcpy(&storage.chip_arr[*way][*chnl].data_blk[*sb][start_blk].page_arr[start_idx], valid_list, sizeof(Page) * dat);

                remain_valid -= dat;
                v += dat;
            }
            else {
                if(remain_valid > PAGE_PER_BLOCK) {
                    memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][start_blk+j].page_arr, &valid_list[v], sizeof(Page) * PAGE_PER_BLOCK);

                    remain_valid -= PAGE_PER_BLOCK;
                    v += PAGE_PER_BLOCK;
                }
                else {
                    memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][start_blk+j].page_arr, &valid_list[v], sizeof(Page) * remain_valid);

                    remain_valid = 0;
                }
            }

            j++;
        }
        
        // Manage ppa_flag
        int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
        memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SBLK);

        memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK * evict_blk], sizeof(int) * (PAGE_PER_BLOCK * (SBLK_SZ - evict_blk)));
        free(global_symb.ppa_flag[superblk]);
        
        global_symb.ppa_flag[superblk] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
        memset(global_symb.ppa_flag[superblk], 0, sizeof(int) * PAGE_PER_SBLK);
        memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SBLK);

        symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, (PAGE_PER_BLOCK * start_blk) + start_idx, num_valid, *sb);
        
        disk_write_cnt += num_valid;

        //free(tmp);
        free(tmp_ppa_flag);
        free(valid_list);
        free(candi_list);
    }
}

void bloom_gc(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, Page* evict_list, int num_evict, int* sb) {
	Page *candi_list;
	int num_candi = 0;
	Page *valid_list, prev_page;
	uint32_t hashkey;
	int num_valid = 0;

	if(!evict_list) {
		int superblk = (*chip) * SBLK_PER_CHIP + (*sb);

		candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
		valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

		memset(candi_list, 0, sizeof(Page) * PAGE_PER_BLOCK);
		memset(valid_list, 0, sizeof(Page) * PAGE_PER_BLOCK);

		// Copy valid pages in evict block
		for(int i=PAGE_PER_BLOCK-1; i>=0; i--) {
			bool is_valid = true;

			for(int b=SBLK_SZ-1; b>0; b--) {
				for(int p=PAGE_PER_BLOCK-1; p>=0; p--) {
					if(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i].oob == 
							storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
						is_valid = false;
						break;
					}
				}

				if(!is_valid) {
					break;
				}
			}

			if(is_valid == true) {
				candi_list[num_candi++] = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr[i];
			}
		}

		qsort(candi_list, num_candi, sizeof(Page), lba_compare);

		// Remove duplicate
		prev_page = valid_list[0] = candi_list[0];
		num_valid++;
		for(int i=1; i<num_candi; i++) {
			if(prev_page.oob != candi_list[i].oob) {
				valid_list[num_valid++] = candi_list[i];
			}

			prev_page = candi_list[i];
		}

		// Free exist block
		free(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr);
		evict_blk_cnt++;

		// Re-allocate and Initialize new block
		Page* tmp_page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

		storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

		memset(storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr, 0, sizeof(Page) * PAGE_PER_BLOCK);
		storage.chip_arr[*way][*chnl].empty[*sb] = (3 * PAGE_PER_BLOCK);

		int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
		memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SBLK);

		memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK], sizeof(int) * PAGE_PER_BLOCK * 3);

		free(global_symb.ppa_flag[superblk]);

		global_symb.ppa_flag[superblk] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
		memset(global_symb.ppa_flag[superblk], 0, sizeof(int) * PAGE_PER_SBLK);
		memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SBLK);

		tmp_page_arr = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;
		for(int b=0; b<SBLK_SZ-1; b++) {
			storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
				storage.chip_arr[*way][*chnl].data_blk[*sb][b+1].page_arr;
		}
		storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr = tmp_page_arr;

#if DEBUG
        for(int i=0; i<num_valid; i++) {
            printf("valid_list[%d]: %u\n", i, valid_list[i].oob);
        }
#endif

		// Flush valid pages
		memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr, valid_list, sizeof(Page) * num_valid);
		symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, PAGE_PER_BLOCK * 3, num_valid, *sb);

		disk_write_cnt += num_valid;

		//free(tmp_page_arr);
		free(tmp_ppa_flag);
		free(candi_list);
	}
	else { // During RB, GC triggererd
		int superblk = (*chip) * SBLK_PER_CHIP + (*sb);
		int evict_blk = -1;

		candi_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SBLK * 2);
		valid_list = (Page*)malloc(sizeof(Page) * PAGE_PER_SBLK);

		memset(candi_list, 0, sizeof(Page) * PAGE_PER_SBLK * 2);

		// Copy valid pages in evict block
		memcpy(candi_list, evict_list, sizeof(Page) * num_evict);
		num_candi += num_evict;

		do {
			num_valid = 0;
			memset(valid_list, 0, sizeof(Page) * PAGE_PER_SBLK);

			evict_blk++;

			// in evict block
			for(int p=0; p<PAGE_PER_BLOCK; p++) {
				bool is_valid = true;

				if(global_symb.ppa_flag[superblk][p] == 1) {
					for(int i=0; i<num_candi; i++) {
						if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[p].oob == candi_list[i].oob) {
							is_valid = false;
							break;
						}
					}
				}
				else {
					if(p != 0) {
						if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[p-1].oob + 1 == storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[p].oob) {
							for(int i=0; i<num_candi; i++) {
								if(storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[p].oob == candi_list[i].oob) {
									is_valid = false;
									break;
								}
							}
						}
						else {
							is_valid = false;
						}
					}
				}

				if(is_valid == true) {
					candi_list[num_candi++] = storage.chip_arr[*way][*chnl].data_blk[*sb][evict_blk].page_arr[p];
				}
			}

			qsort(candi_list, num_candi, sizeof(Page), lba_compare);

			prev_page = valid_list[0] = candi_list[0];
			num_valid++;
			for(int i=1; i<num_candi; i++) {
				if(prev_page.oob != candi_list[i].oob) {
					valid_list[num_valid++] = candi_list[i];
				}

				prev_page = candi_list[i];
			}

		} while((((SBLK_SZ - 1 - evict_blk) * PAGE_PER_BLOCK) + num_valid > PAGE_PER_SBLK) && (evict_blk + 1 < SBLK_SZ));

		// Free exist block
		evict_blk++;
		evict_blk_cnt += evict_blk;
		int exist_blk_num = SBLK_SZ - evict_blk;

		for(int b=0; b<evict_blk; b++) {
			free(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr);

			storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);
			memset(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr, 0, sizeof(Page) * PAGE_PER_BLOCK);
		}

		if(evict_blk == SBLK_SZ) {
			storage.chip_arr[*way][*chnl].empty[*sb] = 0;
		}
		else {
			storage.chip_arr[*way][*chnl].empty[*sb] -= (PAGE_PER_BLOCK * evict_blk);
		}

		// Re-allocate and Initialize new block
		Page* tmp_page_arr = (Page*)malloc(sizeof(Page) * PAGE_PER_BLOCK);

		int* tmp_ppa_flag = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
		memset(tmp_ppa_flag, 0, sizeof(int) * PAGE_PER_SBLK);

		memcpy(tmp_ppa_flag, &global_symb.ppa_flag[superblk][PAGE_PER_BLOCK * evict_blk], \
				sizeof(int) * PAGE_PER_BLOCK * exist_blk_num);

		free(global_symb.ppa_flag[superblk]);

		global_symb.ppa_flag[superblk] = (int*)malloc(sizeof(int) * PAGE_PER_SBLK);
		memset(global_symb.ppa_flag[superblk], 0, sizeof(int) * PAGE_PER_SBLK);
		memcpy(global_symb.ppa_flag[superblk], tmp_ppa_flag, sizeof(int) * PAGE_PER_SBLK);

		for(int i=0; i<evict_blk; i++) {
			tmp_page_arr = storage.chip_arr[*way][*chnl].data_blk[*sb][0].page_arr;
			for(int b=0; b<SBLK_SZ-1; b++) {
				storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr =
					storage.chip_arr[*way][*chnl].data_blk[*sb][b + 1].page_arr;
			}
			storage.chip_arr[*way][*chnl].data_blk[*sb][3].page_arr = tmp_page_arr;
		}

		// Flush valid pages
		int remain_valid = num_valid;
		int j = 0; int v = 0;
		int new_ppn = storage.chip_arr[*way][*chnl].empty[*sb] % PAGE_PER_BLOCK;
		int append_blk = storage.chip_arr[*way][*chnl].empty[*sb] / PAGE_PER_BLOCK;

		while(remain_valid > 0) {
			if(j == 0) {
				memcpy(&storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk].page_arr[new_ppn], valid_list, sizeof(Page) * (PAGE_PER_BLOCK - new_ppn));

				remain_valid -= (PAGE_PER_BLOCK - new_ppn);
				v += (PAGE_PER_BLOCK - new_ppn);

				disk_write_cnt += (PAGE_PER_BLOCK - new_ppn);
			}
			else {
				if(remain_valid > PAGE_PER_BLOCK) {
					memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk + j].page_arr, &valid_list[v], sizeof(Page) * PAGE_PER_BLOCK);

					remain_valid -= PAGE_PER_BLOCK;
					v += PAGE_PER_BLOCK;

					disk_write_cnt += PAGE_PER_BLOCK;
				}
				else {
					if(remain_valid < 0) {
						break;
					}

					memcpy(storage.chip_arr[*way][*chnl].data_blk[*sb][append_blk + j].page_arr, &valid_list[v], sizeof(Page) * remain_valid);

					remain_valid = 0;

					disk_write_cnt += remain_valid;
				}
			}

			j++;
		}
		symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, PAGE_PER_BLOCK * append_blk + new_ppn, num_valid, *sb);

		disk_write_cnt += num_valid;

		//free(tmp_page_arr);
		free(tmp_ppa_flag);
		free(candi_list);
	}

	free(valid_list);
}

void bloom_rebloom_revision(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, int* sb) {
    int superblk = (*chip) * SBLK_PER_CHIP + (*sb);
    Page *candi_list, *rebloom_list, prev_page;
    int num_candi = storage.chip_arr[*way][*chnl].empty[*sb];
    int num_rebloom = 0, candi_sz = 0;
    int gidx = 0;

    candi_list = (Page*)malloc(sizeof(Page) * num_candi);
    rebloom_list = (Page*)malloc(sizeof(Page) * num_candi);

    // Copy data whose ppa_flag is 1
    for(int b=0; b<=(*blk); b++) {
        for(int p=0; p<PAGE_PER_BLOCK; p++) {
            gidx = b * PAGE_PER_BLOCK + p;

            if(global_symb.ppa_flag[superblk][gidx] == 1) {
                candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
            }
            else {
                if(p != 0) {
                    if(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p-1].oob + 1 ==
                            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                        candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
                    }
                }
                else if(b != 0) {
                    if(storage.chip_arr[*way][*chnl].data_blk[*sb][b-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 ==
                            storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                        candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
                    }
                }
                else {
                    candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
                }
            }

            if(gidx == num_candi - 1) {
                break;
            }
        }
    }
    
    qsort(candi_list, candi_sz, sizeof(Page), lba_compare);

    // Make rebloomed-list
    prev_page = candi_list[0];
    for(int i=1; i<candi_sz; i++) {
        if(prev_page.oob == candi_list[i].oob) {
            if(num_rebloom == 0) {
                rebloom_list[num_rebloom++] = candi_list[i];
            }
            else if(rebloom_list[num_rebloom-1].oob == candi_list[i].oob) {
                rebloom_list[num_rebloom-1] = candi_list[i];
            }
            else {
                rebloom_list[num_rebloom++] = candi_list[i];
            }
        }
        else if(prev_page.oob + 1 == candi_list[i].oob) {
            if(num_rebloom == 0) {
                rebloom_list[num_rebloom++] = prev_page;
            }
            else if(rebloom_list[num_rebloom-1].oob == prev_page.oob) {
                rebloom_list[num_rebloom-1] = prev_page;
            }
            else {
                rebloom_list[num_rebloom++] = prev_page;
            }

            rebloom_list[num_rebloom++] = candi_list[i];
        }

        prev_page = candi_list[i];
    }

    // Set invalid bit
    for(int i=0; i<num_rebloom; i++) {
        for(int b=0; b<=(*blk); b++) {
            for(int p=0; p<PAGE_PER_BLOCK; p++) {
                gidx = b * PAGE_PER_BLOCK + p;

                if(rebloom_list[i].oob == storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
                    global_symb.ppa_flag[superblk][gidx] = 0;
                }
                
                if(gidx == num_candi - 1) {
                    break;
                }
            }
        }
    }

    if(num_candi + num_rebloom >= PAGE_PER_SBLK) {
        bloom_gc_revision(pbn, chip, way, chnl, blk, ppn, rebloom_list, num_rebloom, sb);

        free(rebloom_list);
        free(candi_list);
        return;
    }

    // Append rebloomed data
    for(int i=0, p=num_candi%PAGE_PER_BLOCK, b=(*blk); i<num_rebloom; i++, p++) {
        if(p == PAGE_PER_BLOCK) {
            b++;
            p = 0;
        }

        storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p] = rebloom_list[i];
        disk_write_cnt++;
    }

    // resymbolize
    symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, num_rebloom, *sb);
    free(rebloom_list);
    free(candi_list);
}

void bloom_rebloom(uint32_t* pbn, int* chip, int* way, int* chnl, int* blk, int* ppn, int* sb) {
	Page *candi_list, *evict_list, prev_page;
	uint32_t hashkey;
	int superblk = (*chip) * SBLK_PER_CHIP + (*sb);
	
    int num_candi = 0, num_evict = 0, candi_sz = 0;
	num_candi += storage.chip_arr[*way][*chnl].empty[*sb];

	candi_list = (Page*)malloc(sizeof(Page) * num_candi);
	evict_list = (Page*)malloc(sizeof(Page) * num_candi);

	// Copy valid pages and Sort with LBA order
	for(int b=0; b<=(*blk); b++) {
		for(int p=0; p<PAGE_PER_BLOCK; p++) {
			int gidx = b * PAGE_PER_BLOCK + p;

			if(global_symb.ppa_flag[superblk][gidx] == 1) {
				candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
			}
			else {
				if(p != 0) {
					if(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p-1].oob + 1 ==
							storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob) {
						candi_list[candi_sz++] = storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p];
					}
				}
			}

			if(gidx == num_candi - 1) {
				break;
			}
		}
	}

	qsort(candi_list, candi_sz, sizeof(Page), lba_compare);

	// Elect rebloomed-pages
	prev_page = candi_list[0];
	for(int i=1; i<candi_sz; i++) {
		if(prev_page.oob == candi_list[i].oob) {
			if(num_evict == 0) {
				evict_list[num_evict++] = candi_list[i];
			}
			else {
				if(evict_list[num_evict - 1].oob == candi_list[i].oob) {
					evict_list[num_evict - 1] = candi_list[i];
				}
				else {
					evict_list[num_evict++] = candi_list[i];
				}
			}
		}
		else if(prev_page.oob + 1 == candi_list[i].oob) {
			if(num_evict == 0) {
				evict_list[num_evict++] = prev_page;
				evict_list[num_evict++] = candi_list[i];
			}
			else {
				if(evict_list[num_evict - 1].oob == prev_page.oob) {
					evict_list[num_evict - 1] = prev_page;
					evict_list[num_evict++] = candi_list[i];
				}
				else {
					evict_list[num_evict++] = prev_page;
					evict_list[num_evict++] = candi_list[i];
				}
			}
		}

		prev_page = candi_list[i];
	}

	// Set invalid bit with elected rebloomed-pages
	for(int i=0; i<num_evict; i++) {
		for(int b=0; b<(*blk); b++) {
			for(int p=0; p<PAGE_PER_BLOCK; p++) {
				int gidx = b * PAGE_PER_BLOCK + p;

				if(storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p].oob == evict_list[i].oob) {
					global_symb.ppa_flag[superblk][gidx] = 0;
				}

				if(gidx == num_candi - 1) {
					break;
				}
			}
		}
	}

	if(num_candi + num_evict >= PAGE_PER_SBLK) {
		bloom_gc(pbn, chip, way, chnl, blk, ppn, evict_list, num_evict, sb);

		free(evict_list);
		free(candi_list);
		return;
	}

#if DEBUG
	printf("\nIN RB, evict list\n");
	for(int i=0; i<num_evict; i++) {
		printf("evict_list[%d]: %u\n", i, evict_list[i].oob);
	} printf("\n");
	fflush(stdout);
#endif

	// Append rebloomed pages
	int remain_evict = num_evict;
	for(int i=0, p=num_candi%PAGE_PER_BLOCK, b=(*blk); i<num_evict; i++, p++) {
		if(p == PAGE_PER_BLOCK) {
			b++;
			p = 0;
		}

		storage.chip_arr[*way][*chnl].data_blk[*sb][b].page_arr[p] = evict_list[i];
		disk_write_cnt++;
	}

	symbol_resymbolize(*pbn, *chip, *way, *chnl, *blk, num_candi, remain_evict, *sb);
	free(evict_list);
	free(candi_list);
}

void debugging(int way, int chnl, int blk, int pbn) {
	int sum_num_bf = 0;
	uint64_t bfsum = 0, targetsize;

	printf("ppn\tppa\twrite_lba\tppa_flag(BF exists)\tlba_flag(lba exists)\tlba_list\n");

	int chip = way * CHANNEL + chnl;

	for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
		int superblk = chip * SBLK_PER_CHIP + sb;
		sum_num_bf = 0;

		if(superblk != sb) {
			continue;
		}
		printf("[SUPERBLOCK %d in chip %d]\n", sb, way*CHANNEL+chnl);

		for(int b=0; b<SBLK_SZ; b++) {
			printf("[Block %d]\n", b);
			printf("idx\twrite_lba\tppa_flag(BF)\tlba_flag\tlbalist\t\tnew_bf_idx\n");

			for(int p=0; p<PAGE_PER_BLOCK; p++) {
				int pa = b * PAGE_PER_BLOCK + p;

				printf("%d\t\t%u\t\t\t", p, storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[p].oob);
				printf("%d\t\t\t\t%u\t\t%d\t\t\t", global_symb.ppa_flag[superblk][pa], \
						global_symb.lba_flag[superblk][pa], pa);

				if(global_symb.ppa_flag[superblk][pa] == 1 && pa != 0) {
					sum_num_bf++;
					printf("%d\n", sum_num_bf);
				}
				else {
					printf("\n");
				}
			} printf("\n");
		} printf("\n");

		printf("Superblock %d (# of valid BFs: %d)\n", sb, sblk_man->num_bf[superblk]);
	} printf("\n");
}

void bloom_write(uint32_t lba, uint32_t value, char* pttr) {
	uint32_t pbn, hashkey;
	uint32_t superblk;
	int chip, way, chnl, blk, ppn;
	int lba_start, new_bf_idx;
	Page prev_page;
	int sb;

	superblk = lba / mask;
	sb = superblk % (BLOCK_PER_CHIP / SBLK_SZ);

	chip = superblk / SBLK_PER_CHIP;
	way = chip / CHANNEL;
	chnl = chip % CHANNEL;

	ppn = storage.chip_arr[way][chnl].empty[sb]; // TODO: global ppn
	blk = ppn / PAGE_PER_BLOCK; // in one superblock
	pbn = sb * SBLK_SZ + blk; // in one chip

    // GC triggered
	if(ppn >= PAGE_PER_SBLK) {
#if DEBUG
		printf("\nBefore GC in block %u of superblock %u (# of valid BFs: %d)\n", blk-1, \
				superblk, sblk_man->num_bf[superblk]);
		debugging(way, chnl, blk, pbn);

		printf("Before GC ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);
#endif
		bloom_gc_revision(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0, &sb);

		ppn = storage.chip_arr[way][chnl].empty[sb];
		blk = ppn / PAGE_PER_BLOCK;
		pbn = sb * SBLK_SZ + blk;

#if DEBUG
		printf("After GC ppn: %d\tblk: %d\tpbn: %d\n\n", ppn, blk, pbn); fflush(stdout);

		printf("\nAfter GC in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
				superblk, sblk_man->num_bf[superblk]);
		debugging(way, chnl, blk, pbn);
#endif
	}

	// RB triggered (# of valid BFs in one superblock exceeds 50%)
	if(sblk_man->num_bf[superblk] >= REBLOOM) {
#if DEBUG
		printf("\nBefore RB in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
				superblk, sblk_man->num_bf[superblk]);
		debugging(way, chnl, blk, pbn);

		printf("Before RB ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);
#endif

		bloom_rebloom_revision(&pbn, &chip, &way, &chnl, &blk, &ppn, &sb);

		ppn = storage.chip_arr[way][chnl].empty[sb];
		blk = ppn / PAGE_PER_BLOCK;
		pbn = sb * SBLK_SZ + blk;

#if DEBUG
		printf("After RB ppn: %d\tblk: %d\tpbn: %d\n", ppn, blk, pbn); fflush(stdout);

		printf("\nAfter RB in block %u of superblock %u (# of valid BFs: %d)\n", blk, \
				superblk, sblk_man->num_bf[superblk]);
		debugging(way, chnl, blk, pbn);
#endif
	}

	int global_lba_start = ppn;
	lba_start = global_lba_start % PAGE_PER_BLOCK;
	new_bf_idx = sblk_man->num_bf[superblk] + 1;
	bool create_new_bf = false;

	storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lba_start].oob = lba;
	storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lba_start].page = value;
	write_cnt++;
	disk_write_cnt++;

	global_symb.lba_flag[superblk][lba % ((int)(PAGE_PER_SBLK * (1 - OP)))] = 1;

	if(global_lba_start != 0) {
		if(lba_start == 0) {
			if(storage.chip_arr[way][chnl].data_blk[sb][blk-1].page_arr[PAGE_PER_BLOCK-1].oob + 1 == lba) {
				global_symb.ppa_flag[superblk][global_lba_start] = 0;
			}
			else {
				create_new_bf = true;
			}
		}
		else {
			if(storage.chip_arr[way][chnl].data_blk[sb][blk].page_arr[lba_start - 1].oob + 1 == lba) {
				global_symb.ppa_flag[superblk][global_lba_start] = 0;
			}
			else {
				create_new_bf = true;
			}
		}
	}
	else {
        global_symb.ppa_flag[superblk][global_lba_start] = 1;
	}

	if(create_new_bf == true) {
		hashkey = hashing_key(lba);
		symbol_set(bf_man->bits_per_pg[new_bf_idx], new_bf_idx, \
				hashkey + new_bf_idx, global_symb.symbol[chip][sb], \
				st_man->sym_start, st_man->sym_bits_pg);

		sblk_man->num_bf[superblk]++;
		global_symb.ppa_flag[superblk][global_lba_start] = 1;
	}

	storage.chip_arr[way][chnl].empty[sb]++;
}

void bloom_read(uint32_t lba) {
    uint32_t pbn;
    register uint32_t hashkey;
    register uint8_t* reg_symbol;
	uint32_t superblk;
	int chip, way, chnl, blk, ppn;
	uint32_t real_lba;
	int locate, new_bf_idx;
	int delta=0, sb;

	superblk = lba / mask;
	sb = superblk % (BLOCK_PER_CHIP / SBLK_SZ);

	chip = superblk / SBLK_PER_CHIP;
	way = chip / CHANNEL;
	chnl = chip % CHANNEL;

	ppn = storage.chip_arr[way][chnl].empty[sb];
	blk = ppn / PAGE_PER_BLOCK;
	pbn = sb * SBLK_SZ + blk;

	locate = lba % ((int)(PAGE_PER_SBLK * (1 - OP)));

    new_bf_idx = sblk_man->num_bf[superblk];

    reg_symbol = global_symb.symbol[chip][sb];

	if(locate != 0) {
		for(int i=locate-1; i>=0; i--) {
			if(global_symb.lba_flag[superblk][i] == 0) {
				break;
			}
			else {
				delta++;
			}
		}
	}

	real_lba = lba - delta;
	hashkey = hashing_key(real_lba);

    int idx = (ppn % PAGE_PER_BLOCK) - 1;
    int pa;
    int seq_cnt = 0;

	for(int b=blk; b>=0; b--) {
		for(; idx>=0; idx--) {
			if(b == 0 && idx == 0) {
                break;
            }

            pa = b * PAGE_PER_BLOCK + idx;

            if(global_symb.ppa_flag[superblk][pa] == 1) {
                if(symbol_check(bf_man->bits_per_pg[new_bf_idx], hashkey + new_bf_idx, \
                            reg_symbol, st_man->sym_start[new_bf_idx], st_man->sym_bits_pg[new_bf_idx]) == true) {
                    if(storage.chip_arr[way][chnl].data_blk[sb][b].page_arr[idx].oob == real_lba) {
                        found_cnt++;
                        return;
                    }
                    else {
                        notfound_cnt++;
                        seq_cnt = 0;
                    }
                }
                else {
                    false_cnt++;
                    seq_cnt = 0;
                }

                new_bf_idx--;
            }
            else {
                seq_cnt++;
            }

            read_loop++;
		}

        idx = PAGE_PER_BLOCK - 1;
	}

    // Case of page offset 0
    if(storage.chip_arr[way][chnl].data_blk[sb][0].page_arr[0].oob == real_lba) {
        found_cnt++;
    }
    else {
        printf("\nWRITE IS STOPPED\n");
        debugging(way, chnl, blk, pbn);
        printf("This should not happen !!\n");
        printf("Until now, found cnt: %d (LBA %u delta %d)\n", found_cnt, lba, delta);
        uint32_t startlba = lba;
        for(int i=4; i>0; i--) {
            printf("lbalist %u\tlba_flag %d\n", startlba - i, global_symb.lba_flag[superblk][(startlba - i) % ((int)(PAGE_PER_SBLK * (1 - OP)))]);
        }
        for(int i=0; i<5; i++) {
            printf("lbalist %u\tlba_flag %d\n", startlba + i, global_symb.lba_flag[superblk][(startlba + i) % ((int)(PAGE_PER_SBLK * (1 - OP)))]);
        }
        printf("\n");
        fflush(stdout);
        exit(1);
    }
}

#if BUFFERED
void bloom_write(uint32_t lba, uint32_t value, char* pttr) {
	uint32_t pbn, hashkey;
	int chip, way, chnl, blk, ppn;
	int lba_start, buf_sz, cur_num_bf;
	Page prev_page;

	pbn = lba >> mask; // random write
	chip = pbn / BLOCK_PER_CHIP;
	way = chip / CHANNEL;
	chnl = chip % CHANNEL;
	blk = pbn % BLOCK_PER_CHIP;

	ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

	if(ppn >= BLOCK_PPN) {
		bloom_gc(&pbn, &chip, &way, &chnl, &blk, &ppn, NULL, 0);
	}

	// The number of valid BFs in one superblock exceeds 50%
	if(sblk_man->num_bf[pbn] >= REBLOOM) {
		bloom_rebloom(&pbn, &chip, &way, &chnl, &blk, &ppn);
		ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;
	}

	buf_sz = block_buffer[chip][blk].buf_sz;

	// Buffering LBA
	if(buf_sz < BUF_SZ) {
		block_buffer[chip][blk].lba[buf_sz] = lba;
		block_buffer[chip][blk].value[buf_sz] = value;
		block_buffer[chip][blk].buf_sz++; buf_sz++;
	}

	// Flush write buffer
	if(buf_sz == BUF_SZ) {
		lba_start = ppn * BUF_SZ;
		cur_num_bf = sblk_man->num_bf[pbn];

		if(lba_start != 0) {
			prev_page = storage.chip_arr[way][chnl].data_blk[blk].page_arr[lba_start - 1];
		}

		// TODO: how to synchronize BF's start? 0 or 1?
		for(int i=0; i<BUF_SZ; i++) {
			if(lba_start + i != 0) {
				if(prev_page.oob + 1 == block_buffer[chip][blk].lba[i]) {
					global_symb.ppa_flag[pbn][lba_start + i] = 0;
				}
				else {
					hashkey = hashing_key(block_buffer[chip][blk].lba[i]);
					symbol_set(bf_man->bits_per_pg[cur_num_bf + i], cur_num_bf + i, \
							hashkey + cur_num_bf + i, global_symb.symbol[chip][blk], \
							st_man->sym_start, st_man->sym_bits_pg);

					sblk_man->num_bf[pbn]++;
					global_symb.ppa_flag[pbn][lba_start + i] = 1;
				}
			}

			storage.chip_arr[way][chnl].data_blk[blk].page_arr[lba_start + i].page = block_buffer[chip][blk].value[i];
			storage.chip_arr[way][chnl].data_blk[blk].page_arr[lba_start + i].oob = block_buffer[chip][blk].lba[i];

			global_symb.lba_flag[pbn][block_buffer[chip][blk].lba[i] % PAGE_PER_BLOCK] = 1;
			prev_page.oob = block_buffer[chip][blk].lba[i];
			disk_write_cnt++;
		}

		storage.chip_arr[way][chnl].data_blk[blk].empty++;
		buf_sz = block_buffer[chip][blk].buf_sz = 0;

		write_cnt += BUF_SZ;
	}
}

void bloom_read(uint32_t lba) {
	uint32_t pbn;
    register uint32_t hashkey;
    register uint8_t* reg_symbol;
	int chip, way, chnl, blk, ppn;
	struct timeval strt, end;

	bool rb_or_gc;
	int cur_key = lba % PAGE_PER_BLOCK;
	uint32_t real_lba;
	int delta=0, new_bf_idx;

	pbn = lba >> mask;
	chip = pbn / BLOCK_PER_CHIP;
	way = chip / CHANNEL;
	chnl = chip % CHANNEL;
	blk = pbn % BLOCK_PER_CHIP;

	ppn = storage.chip_arr[way][chnl].data_blk[blk].empty;

	new_bf_idx = sblk_man->num_bf[pbn];
	rb_or_gc = (new_bf_idx == (ppn * BUF_SZ)) ? false : true;

    reg_symbol = global_symb.symbol[chip][blk];

	if(rb_or_gc == true) {
		// Get sequential delta first
		if(cur_key != 0) {
			for(int idx=cur_key-1; idx>=0; idx--) {
				if(global_symb.lba_flag[pbn][idx] == 0) {
					break;
				}
				else {
					delta++;
				}
			}
		}

		real_lba = lba - delta;
		hashkey = hashing_key(lba - delta);

		int seq_cnt = 0;
		for(int idx=ppn*BUF_SZ-1, j=new_bf_idx; idx>0; idx--) {
			if(global_symb.ppa_flag[pbn][idx] == 1) {
				if(symbol_check(bf_man->bits_per_pg[j], hashkey + j, \
							reg_symbol, st_man->sym_start[j], st_man->sym_bits_pg[j]) == true) {
					if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == real_lba) {
						found_cnt++;
						return;
					}
					else {
						notfound_cnt++;
					}
				}
				else {
					false_cnt++;
				}

				j--;
			}

			read_loop++;
		}
	}
	else {
		hashkey = hashing_key(lba);

		for(int idx=ppn*BUF_SZ-1; idx>0; idx--) {
			if(symbol_check(bf_man->bits_per_pg[idx], hashkey + idx, \
						reg_symbol, st_man->sym_start[idx], st_man->sym_bits_pg[idx]) == true) {
				if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[idx].oob == lba) {
					found_cnt++;
					return;
				}
				else {
					notfound_cnt++;
				}
			}
			else {
				false_cnt++;
			}

			read_loop++;
		}
	}

	// Case of page offset 0
	if(storage.chip_arr[way][chnl].data_blk[blk].page_arr[0+delta].oob != lba) {
		printf("This should not happen !!\n");
		exit(1);
	}
	else {
		found_cnt++;
	}
}
#endif

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

		st_man->sym_start[p] = sum;
		sum += symbol_bit;
	}
	
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

void print_stats(char* w_type, char* r_type) {
	uint64_t sum = 0;
	uint64_t targetsize = 0;

	printf("\nTEST DATASET: write %d read %d\n", DATA_SET, DATA_SET);
	printf("TEST TYPE: %s write, %s read\n", w_type, r_type);
	
    printf("\n### BENCHMARK RESULTS - READ ###\n");
	printf("Total read requests: %d\n", read_cnt);
	printf("Total found num: %d\n", found_cnt);
	printf("Total not-found num: %d\n", notfound_cnt);
	printf("Total false num: %ld\n", false_cnt);
	printf("RAF: %.2f\n", (double)(found_cnt + notfound_cnt) / read_cnt);
	
    printf("\n### BENCHMARK RESULTS - WRITE ###\n");
	printf("Total write requests: %d\n", write_cnt);
	printf("Total number of disk write: %d\n", disk_write_cnt);
	printf("WAF: %.2f\n", (double)disk_write_cnt / write_cnt);
	printf("Total erase count: %d\n", evict_blk_cnt);

	printf("\n### BLOOMFILTER INFO ###\n");
	printf("Sum of BF bits in 1 Superblock: ");
	for(int p=0; p<PAGE_PER_SBLK; p++) {
		sum += bf_man->bits_per_pg[p];
	}
	targetsize = sum / 8;
	if(sum % 8) {
		targetsize++;
	}
	printf("%lu bits, %lu bytes (%.2lf%% of PFTL: BIT UNIT)\n", sum, targetsize, (double)sum/BIT_PFTL*100);

	uint64_t bytes_per_chip = 0;
	printf("Sum of BF bytes in 1 Chip: ");
	bytes_per_chip += (targetsize * SBLK_PER_CHIP);
	printf("%lu bytes (%.2lf%% of PFTL: BYTE UNIT)\n", bytes_per_chip, (double)bytes_per_chip/(BYTE_PFTL*SBLK_PER_CHIP)*100);

	uint64_t total_bytes = 0;
	printf("Total BF bytes: ");
	total_bytes += (bytes_per_chip * CHIP);
	printf("%lu bytes (%.2lf%% of PFTL: BYTE UNIT)\n", total_bytes, (double)total_bytes/(BYTE_PFTL*SBLK_PER_CHIP*CHIP)*100);

	printf("\n### SYMBOL TABLE INFO ###\n");
	targetsize = st_man->sym_bits_super_blk / 8;
	if(st_man->sym_bits_super_blk % 8) {
		targetsize++;
	}
	printf("Sum of SYMB bits in 1 Superblock: %lu bits, %lu bytes", st_man->sym_bits_super_blk, targetsize);
	printf(" (%.2lf%% of PFTL: BIT UNIT)\n", (double)st_man->sym_bits_super_blk/BIT_PFTL*100);

	targetsize = st_man->sym_bits_chip / 8;
	if(st_man->sym_bits_chip % 8) {
		targetsize++;
	}
	printf("Sum of SYMB bytes in 1 Chip: %lu bytes ", targetsize);
    printf("(%.2lf%% of PFTL: BYTE UNIT)\n", (double)targetsize/(BYTE_PFTL*SBLK_PER_CHIP)*100);

	targetsize = st_man->sym_bits_total / 8;
	if(st_man->sym_bits_total % 8) {
		targetsize++;
	}
	printf("Total SYMB bytes: %lu bytes ", targetsize);
    printf("(%.2lf%% of PFTL: BYTE UNIT)\n", (double)targetsize/(BYTE_PFTL*SBLK_PER_CHIP*CHIP)*100);

	sum = 0;
	printf("\n### AFTER REBLOOMING ###\n");
    printf("Average number of BF: ");
    for(int b=0; b<TOTAL_SBLK; b++) {
        sum += sblk_man->num_bf[b];
    }
    printf("%ld\n", sum/TOTAL_SBLK);

	for(int p=0; p<=sblk_man->num_bf[0]; p++) {
		sum += st_man->sym_bits_pg[p];
	}
	sum += (2 * PAGE_PER_SBLK);
	targetsize = sum / 8;
	if(sum % 8) {
		targetsize++;
	}
	printf("Sum of REBLOOMED bits in 1 Superblock: %lu bits, %lu bytes ", sum, targetsize);
	printf("(%.2lf%% of PFTL: BIT UNIT)\n", (double)sum/BIT_PFTL*100);

    sum = 0;

	sum = 0;
	for(int b=0; b<TOTAL_SBLK; b++) {
		for(int p=0; p<=sblk_man->num_bf[b]; p++) {
			sum += st_man->sym_bits_pg[p];
		}

		sum += (2 * PAGE_PER_SBLK);
	}
    targetsize = sum / 8;
    if(sum % 8) {
        targetsize++;
    }
    printf("Total REBLOOMED bytes: %lu bytes ", targetsize);
    printf("(%.2lf%% of PFTL: BYTE UNIT)\n", (double)targetsize/(BYTE_PFTL*SBLK_PER_CHIP*CHIP)*100);

    printf("\n### TIME RECORDS - READ ###\n");
    printf("Total read time: %.f (us)\n", r_time);
	printf("Average read time: %f (us) *****\n", r_time/DATA_SET);
	printf("Total read loop count: %d\n", read_loop);
	printf("Total read check time: %.f (us)\n", read_check);
	printf("Read check time per req: %f (us)\n", read_check/DATA_SET);
	printf("Average read check time: %f (us)\n", read_check/read_loop);
    
    printf("\n### TIME RECORDS - WRITE ###\n");
	printf("Total write time: %.f (us)\n", w_time);
	printf("Average write time: %f (us)\n", w_time/DATA_SET);

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

	evict_blk_cnt = read_cnt = write_cnt = true_cnt = false_cnt = found_cnt = notfound_cnt = 0;

	if(argc != 3) {
		printf("Usage: ./simulationFTL W-type R-type\n(Type must be either SEQ or RAND)\n");
		exit(1);
	}

	bloom_init();
	symbol_init();

	// Allocate LBA array
	read_arr = (uint32_t*)malloc(sizeof(uint32_t) * DATA_SET);
	write_arr = (uint32_t*)malloc(sizeof(uint32_t) * DATA_SET);

	printf("\nTEST START !!\n");
	srand((unsigned int)time(NULL));

	// Generate LBA list
	// Options: W_UNIQUE | R_UNIQUE | R_RD_ORDER | R_RD_GEN
	op = R_RD_GEN;
	//op = W_UNIQUE | R_UNIQUE;
	make_test_set(write_arr, read_arr, argv[1], argv[2], op, DATA_SET);

	val = 0;
	while(try_wr--) {
		lba = write_arr[val];

		gettimeofday(&strt, NULL);
		bloom_write(lba, val, argv[1]);

		gettimeofday(&end, NULL);
		w_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));

		val++;
	}

	int try = HEAVY_GC;
	while(try--) {
		try_wr = DATA_SET; 
		val = 0;

		memset(write_arr, 0, sizeof(uint32_t) * DATA_SET);
		make_test_set(write_arr, read_arr, "RAND", argv[2], op, DATA_SET);

		while(try_wr--) {
			lba = write_arr[val];

			gettimeofday(&strt, NULL);
			bloom_write(lba, val, argv[1]);
			gettimeofday(&end, NULL);
			w_time += (((end.tv_sec - strt.tv_sec) * 1000000) + (end.tv_usec - strt.tv_usec));

			val++;
		}
	}

#if DEBUG
	printf("\nWRITE DONE\n");
	int chip = 0;

	for(int sb=0; sb<SBLK_PER_CHIP; sb++) {
		int superblk = chip * SBLK_PER_CHIP + sb;
		int sum_num_bf = 0;

		printf("[SUPERBLOCK %d in chip %d]\n", sb, 0);

		for(int b=0; b<SBLK_SZ; b++) {
			printf("[Block %d]\n", b);
			printf("idx\twrite_lba\tppa_flag(BF)\tlba_flag\tlbalist\t\tnew_bf_idx\n");

			for(int p=0; p<PAGE_PER_BLOCK; p++) {
				int pa = b * PAGE_PER_BLOCK + p;

				printf("%d\t\t%u\t\t\t", p, storage.chip_arr[0][0].data_blk[sb][b].page_arr[p].oob);
				printf("%d\t\t\t\t%u\t\t%u\t\t\t", global_symb.ppa_flag[superblk][pa], \
						global_symb.lba_flag[superblk][pa], pa);

				if(global_symb.ppa_flag[superblk][pa] == 1) {
					sum_num_bf++;
					printf("%d\n", sum_num_bf);
				}
				else {
					printf("\n");
				}
			} printf("\n");
		} printf("\n");

		printf("Superblock %d (# of valid BFs: %d)\n", sb, sblk_man->num_bf[superblk]);
	} printf("\n");
    fflush(stdout);
#endif

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
