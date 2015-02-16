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
#include <list>
#include "seq.h"
#include "klib/khash.h"
#include "gz.h"
#define BUF_STACK (1 << 22)

// twobit, complement-neutral key
#define get_b2cn_key(t, dna, rev) ({\
        t = (dna & KEYNT_STRAND) ? ((rev) ^ 0xaaaaaaaaaaaaaaaa) : dna;\
        (((t >> 1) & ~HALF_KEYNT_MASK) | (t & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;\
})

#define __append_next_b2(b, s, b2pos, dna, rc) ({\
    b = b2pos;\
    ASSERT(b < s ## _l, return -1u, ":%u", b);\
    b = (s[b>>2] >> ((b & 3) << 1)) & 3;\
    rc = ((b ^ 2) << KEYNT_TOP) | (rc >> 2);\
    ((dna << 2) & KEYNT_MASK) | b;\
})

#define __get_next_ndx(ndx, dna, rc) ({\
        ndx = (dna & KEYNT_STRAND) ? dna : rc;\
        ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
})

#define __next_ndx_with_b2(b, s, b2pos, dna, rc) ({\
    dna = __append_next_b2(b, s, b2pos, dna, rc);\
    __get_next_ndx(b, dna, rc);\
})

#define __init_key(i, b, s, b2pos, dna, rc) ({\
    dna = rc = 0; \
    for (i = b2pos + KEY_WIDTH; b2pos != i; ++b2pos){\
        dna = __append_next_b2(b, s, b2pos, dna, rc);\
        /*fputc(b6(b << 1), stderr);*/\
    }\
    /*fputc('\n', stderr);*/\
})

KHASH_MAP_INIT_INT64(UQCT, unsigned)

#define NO_MULTIMAPPER_YET 0
#define NEW_MULTIMAPPER 2
#define _MULTIMAPPER 0xff // 0xff is max for char
#define GSZ_MAX 0xff000000

#define INFERIORITY_BIT (1u<<31)
#define ORIENTATION (1u<<30)

#define N_MASK ((1u << KEY_WIDTH) - 1u)

#define ULL(x) ((unsigned __int128)(x))

#define REF_CHANGE (3u<<30)
#define N_STRETCH (2u<<30)
// ... room for more; KNOWN_POS (dbSNP/cosmic)?
#define RESERVED (1u<<30)

#define NR_MASK 0x3fffffff
#define TYPE_MASK 0xc0000000
#undef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#undef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
typedef struct kcs_t { // Keycounts
    uint32_t b2pos;
    uint32_t infior;
} Kcs;

// Either chromosome boundaries, stretches of N's or regions covered by unique keys.
typedef struct bnd_t {
    uint32_t b2pos;
    uint32_t len;
} Bnd;

typedef struct kct {
    char *hdr;
    uint8_t *seq;
    khash_t(UQCT) *H;
    uint32_t *kcsndx;
    Kcs *kcs;
    std::list<Bnd> bnd;
    unsigned seq_l, kcs_l, bnd_l;
    uint32_t Nmask, hdr_l;
    uint8_t seq_m, hdr_m, kcs_m, kcsndx_m, bnd_m; //XXX: bitfields?
} kct_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
