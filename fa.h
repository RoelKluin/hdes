/*
 *    Description:
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel Kluin,
 */

#ifndef RK_FA_H
#define RK_FA_H
#include "seq.h"
#include "klib/khash.h"
#include "gz.h"
#define BUF_STACK (1 << 22)

// twobit, complement-neutral key
#define get_b2cn_key(t, dna, rev) ({\
        t = (dna & KEYNT_STRAND) ? ((rev) ^ 0xaaaaaaaaaaaaaaaa) : dna;\
        (((t >> 1) & ~HALF_KEYNT_MASK) | (t & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;\
})

KHASH_MAP_INIT_INT(UQCT, unsigned)

typedef struct mmap {
    unsigned* ct;
    unsigned* first_pos;
    uint64_t lastmmap;
    unsigned l, m;
} mmap_t;

typedef struct kct {
    int (*process) (struct seqb2_t*, struct kct*, uint64_t);
    char* x;
    khash_t(UQCT) *H;
    mmap_t* mm;
    uint64_t pos;
    unsigned l;
} kct_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
