#ifndef __SWSHA256_H__
#define __SWSHA256_H__

#include <stdint.h>
#include <string.h>

// calculate SHA-256 hash
// input: uint32_t key
// output: 256 bit result in state[8]
void sha256_calculate(uint32_t state[8], uint32_t key);

#endif
