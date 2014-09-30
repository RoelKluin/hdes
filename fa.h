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
#include "gz.h"
#define BUF_STACK (1 << 22)

#define add_b6N0(t, c, dna, rev) ({\
        t = b6N0(c);\
        dna = (dna << 2) | t;\
        rev = ((uint64_t)t << KEYNT_TOP) | (rev >> 2);\
})

// twobit, complement-neutral key
#define get_b2cn_key(t, dna, rev) ({\
        t = (dna & KEYNT_STRAND) ? ((rev) ^ 0xaaaaaaaaaaaaaaaa) : dna;\
        (((t >> 1) & ~HALF_KEYNT_MASK) | (t & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;\
})

typedef struct index_action {
    const char* search;
    const char* replace;
} index_action_t;

typedef struct kct {
    int (*process) (uint8_t*, struct kct*);
    gzFile out;
    char* x;
    unsigned l;
} kct_t;



int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
