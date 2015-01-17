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

#define N_MASK ((1u << KEY_WIDTH) - 1u)

#define ULL(x) ((unsigned __int128)(x))

typedef struct Kcs {
    unsigned __int128 lastp: 32;
    unsigned __int128 dbits: 5;
    unsigned __int128 res: 3; // reserved (TODO: mode?)
    unsigned __int128 ct1: 28; // 64
    unsigned __int128 ct2: 20;
    unsigned __int128 ct3: 16;
    unsigned __int128 ct4: 24;
} kcs;

typedef struct kct {
    char* x, *hdr;
    uint8_t *seq;
    khash_t(UQCT) *H;
    kcs* kp;
    unsigned* mm, *ho;
    unsigned *at;
    unsigned seq_l, kp_l;
    unsigned at_l; // rep
    unsigned l; // length of char* x;
    uint32_t Nmask, mm_l, hdr_l, ho_l;
    uint8_t mm_m, at_m, seq_m, hdr_m, ho_m, kp_m; //XXX: bitfields?
} kct_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
