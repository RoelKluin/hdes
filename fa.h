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

KHASH_MAP_INIT_INT64(UQCT, unsigned)

#define NO_MULTIMAPPER_YET 0
#define NEW_MULTIMAPPER 2
#define _MULTIMAPPER 0xff // 0xff is max for char
#define GSZ_MAX 0xff000000

#define N_STRETCH 0xfffe
#define REF_CHANGE 0xffff
#define N_MASK ((1u << KEY_WIDTH) - 1u)

#define ULL(x) ((unsigned __int128)(x))

typedef struct Rgn { // Region
    uint32_t pos;
    uint32_t nr:30; // unique iteration count, N-stretch or hdr offset (can decrease bits)
    uint32_t type: 2; // unique pos / N-stretch / Chromo / unique pos already processed
} rgn;

typedef struct Kcs { // Keycounts
    uint32_t lastp;
    uint32_t ct;
} kcs;

typedef struct kct {
    char *hdr;
    uint8_t *seq;
    khash_t(UQCT) *H;
    uint32_t *kpndx;
    kcs *kp;
    rgn *reg;
    unsigned seq_l, kp_l;
    uint32_t Nmask, hdr_l, reg_l;
    uint8_t seq_m, hdr_m, reg_m, kp_m, kpndx_m; //XXX: bitfields?
} kct_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
