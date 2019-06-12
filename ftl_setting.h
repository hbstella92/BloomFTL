#ifndef __BLOOMFTL__
#define __BLOOMFTL__

#include <stdint.h>

#define OP 0.07

// Simulation version
#define DEBUG 0
#define ARM 0
#define REBLOOM (int)((PAGE_PER_SBLK)*0.5)
#define BUFFERED 0
#define HEAVY_GC 0

// FTL
#define CHANNEL 1
#define WAY 1
#define CHIP ((CHANNEL)*(WAY))

// Block
#define BLOCK_PER_CHIP 8
//#define BLOCK_PER_CHIP 8192
#define PAGE_PER_BLOCK 256
#define TOTAL_BLOCK ((CHIP)*(BLOCK_PER_CHIP))

// Superblock
#define SBLK_SZ 4
#define SBLK_PER_CHIP ((BLOCK_PER_CHIP)/(SBLK_SZ))
#define PAGE_PER_SBLK ((PAGE_PER_BLOCK)*(SBLK_SZ))
#define TOTAL_SBLK ((TOTAL_BLOCK)/(SBLK_SZ))
#define TOTAL_PAGE (int)((TOTAL_SBLK)*(int)((PAGE_PER_SBLK)*(1 - OP)))
#define LB_SZ (8*(K))
#define PFTL ((PAGE_PER_SBLK)*32)
#define DATA_SET (TOTAL_PAGE)

#if BUFFERED
#define BUF_SZ 4
#define PP_SZ ((LB_SZ)*(BUF_SZ))
#define PPN_BLK ((PAGE_PER_BLOCK)/(BUF_SZ))
#endif

// LBA generation
#define W_UNIQUE 0x1
#define R_UNIQUE 0x2
#define R_RD_ORDER 0x4
#define R_RD_GEN 0x8

// Unit
#define K (1024)
#define M (1024*K)
#define G (1024*M)
#define T (1024*G)
#define P (1024*T)
#define KEYT uint32_t

#ifndef __GNUG__
typedef enum{false,true} bool;
#endif

#endif
