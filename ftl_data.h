#ifndef __BLOOMDATA__
#define __BLOOMDATA__

#include <stdint.h>

typedef struct {
	uint32_t k;
	uint64_t m;
	uint64_t targetsize;
	int n;
	float p;
	char* body;
	uint64_t start; // start index of BF
} BF;

typedef struct {
	int* num_bf;
} SBlkManager;

typedef struct {
	uint32_t page;
	uint32_t oob;
} Page;

typedef struct {
	Page* page_arr;
} Block;

typedef struct {
	Block** data_blk;
	int* empty;
} Chip;

typedef struct {
	Chip** chip_arr;
} Storage;

typedef struct {
	BF** bfchip_arr;
} GlobalBF;

typedef struct {
	uint64_t* bits_per_pg;
	uint64_t* bytes_arr;
} BFManager;

// Symbol
typedef struct {
	uint8_t*** symbol;

	int** lba_flag;
	int** ppa_flag; // TODO: must change to bit array
} GlobalSymb;

typedef struct {
	uint64_t sym_bytes_total;
	uint64_t dead_bytes_total;

	uint64_t sym_bits_chip; // 1 chip
	uint64_t sym_bits_blk; // 1 block
	uint64_t sym_bits_super_blk; // 1 superblk
	uint64_t* sym_bits_pg;
	uint64_t sym_bits_total;
	uint64_t* new_pidx_start;

	uint64_t* sym_start; // bit unit
} SManager;

#if BUFFERED
typedef struct {
	uint32_t* lba;
	uint32_t* value;
	int buf_sz;
} BlockBuffer; // TODO: change to struct Page
#endif

#endif
