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
#include <stack>
#include <queue>
#include "seq.h"
#include "klib/khash.h"
#include "gz.h"
#define BUF_STACK (1 << 22)

#define _seq_next(b, dna, rc) ({\
    rc = ((b ^ 2) << KEYNT_TOP) | (rc >> 2);\
    ((dna << 2) & KEYNT_MASK) | b;\
})

#define _seq_prev(b, dna, rc) ({\
    rc = ((rc << 2) & KEYNT_MASK) | (b ^ 2);\
    (b << KEYNT_TOP) | (dna >> 2);\
})

#define read_b2(s, p) ((s[(p)>>2] >> (((p) & 3) << 1)) & 3)

#define __rdsq(direction, b, s, b2pos, dna, rc) ({\
    b = b2pos;\
    ASSERT(b < s ## _l, return -1u, ":%u", b);\
    b = read_b2(s, b);\
    _seq_ ## direction (b, dna, rc);\
})

#define _get_ndx(ndx, dna, rc) ({\
    ndx = (dna & KEYNT_STRAND) ? dna : rc;\
    ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
})

#define __rdndx(direction, b, s, b2pos, dna, rc) ({\
    dna = __rdsq(direction, b, s, b2pos, dna, rc);\
    _get_ndx(b, dna, rc);\
})

#define _read_next_ndx(b, s, p, d, r) __rdndx(next, b, s, p, d, r)
#define _read_prev_ndx(b, s, p, d, r) __rdndx(prev, b, s, p - KEY_WIDTH, d, r)

#define _init_key(i, b, s, b2pos, dna, rc) ({\
    dna = rc = 0u; \
    for (i = b2pos + KEY_WIDTH; b2pos != i; ++b2pos){\
        dna = __rdsq(next, b, s, b2pos, dna, rc);\
        /**/if (dbg) { fputc(b6(b << 1), stderr); }/**/\
    }\
    /**/if (dbg) { fputc('\n', stderr); }/**/\
})

KHASH_MAP_INIT_INT64(UQCT, unsigned)

#define UNINITIALIZED (~0u)

#define NO_MULTIMAPPER_YET 0
#define NEW_MULTIMAPPER 2
#define _MULTIMAPPER 0xff // 0xff is max for char
#define GSZ_MAX 0xff000000

#define INFERIORITY_BIT (1u<<31)
#define ORIENTATION (1u<<30)

#define NEEDS_UPDATE 0x8000000

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
    uint16_t infior;
    uint16_t ct;
} Kcs;



typedef struct rng_t {
    uint32_t s, e; //b2pos start and end of range
} Rng;

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
    std::list<Rng> rng;
    unsigned seq_l, kcs_l, bnd_l;
    uint32_t Nmask, hdr_l;
    uint8_t seq_m, hdr_m, kcs_m, kcsndx_m, bnd_m; //XXX: bitfields?
} kct_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
