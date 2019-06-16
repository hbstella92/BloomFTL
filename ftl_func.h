#ifndef __BLOOMFUNC__
#define __BLOOMFUNC__

#include <stdint.h>

#include "ftl_data.h"

void bloom_init();
void bloom_destroy();
void bloom_write(uint32_t, uint32_t, char *);
void bloom_read(uint32_t);
void bloom_gc(uint32_t *, int *, int *, int *, int *, int *, Page *, int, int *);
void bloom_rebloom(uint32_t *, int *, int *, int *, int *, int *, int *);

void print_stats(char *, char *);

static inline uint32_t hashing_key(uint32_t);
int compare(const void *, const void *);
void shuffle(uint32_t *, int);
void make_test_set(uint32_t *, uint32_t *, char *, char *, uint8_t, uint32_t);
void print_byte_as_bit(unsigned char *, size_t);
int lba_compare(const void *, const void *);
void debugging(int, int, int, int);

#endif
