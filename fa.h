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
#include <map>
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

#define _get_ndx(ndx, b, dna, rc) ({\
    b = (dna & KEYNT_STRAND);\
    ndx = b ? dna : rc;\
    ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
})

#define __rdndx(direction, ndx, b, s, b2pos, dna, rc) ({\
    dna = __rdsq(direction, b, s, b2pos, dna, rc);\
    _get_ndx(ndx, b, dna, rc);\
})

#define _read_next_ndx(ndx, b, s, p, d, r) __rdndx(next, ndx, b, s, p, d, r)
#define _read_prev_ndx(ndx, b, s, p, d, r) __rdndx(prev, ndx, b, s, p - KEY_WIDTH, d, r)

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

#define NR_MASK 0x3fffffff
#define TYPE_MASK 0xc0000000
#undef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#undef min
#define min(a,b) ((a) <= (b) ? (a) : (b))

#define UNKNOWN_HDR 0
#define ENSEMBL_HDR 1

#define TEST_BND 1
#ifdef TEST_BND
# define DEBUG_ASSIGN_ENDDNA(end_dna_ref, dna) (end_dna_ref) = (dna)
#define ASSIGN_BD(_bd, _t, _d, _s, _l, _i, _e) \
        _bd = {.t = _t, .dna = _d, .s = _s, .l = _l, .i = _i, .end_dna = _e};
#else
# define DEBUG_ASSIGN_ENDDNA(end_dna_ref, dna) //nothing
#define ASSIGN_BD(_bd, _t, _d, _s, _l, _i, _e) \
        _bd = {.t = _t, .dna = _d, .s = _s, .l = _l, .i = _i};
#endif

typedef packed_struct bnd_t {
    uint64_t t: 8; // type
    uint64_t dna: 56;
    uint32_t s; // start
    uint32_t l; // length
    uint32_t i; // inferiority
#ifdef TEST_BND
    uint64_t end_dna;
#endif
} Bnd;

union Kct {
    packed_struct {
        uint8_t m;
        uint8_t l;
        uint8_t b2[14]; // first 2bits, if exceeding 14*4 2bits, use p below instead.
    } seq;
    packed_struct {
        uint64_t m: 8;
        uint64_t l: 40;
        uint64_t etc: 16;
        uint8_t* b2; // next_b2_0|next_b2_1|next_b2_2|next_b2_3|... 
    } p;
    Bnd rng;
};

packed_struct Walker {
    uint32_t count;
    uint32_t infior: 22; //b2pos start and end of range
    uint32_t tmp_count: 10;
};

struct Hdr {
    uint32_t part[10]; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    uint32_t end_b2pos;
    uint8_t hdr_type;
    uint8_t p_l;
    std::list<uint32_t> bnd; //
};

struct kct_t {
    struct Tid {
        Tid() { id = _buf_init_err(id, 0, id_m = 0x1f); }
        ~Tid() { free(id); }
        int add(char c) {
            _buf_grow_err(id, 1ul, return -ENOMEM);
            id[id_l++] = c;
            return 0;
        }
        char* at(uint32_t l) const { return l < id_l ? id + l : NULL; }
        unsigned l() const { return id_l; }
        unsigned m() const { return id_m; }
        bool operator () (uint32_t a, uint32_t b) const 
        {
            return strcmp(&id[a], &id[b]) == 0;
        }
    private:
        char* id;
        unsigned id_l: 27;
        unsigned id_m: 5;
    } tid;

    uint32_t *kcsndx;
    Kct *kct;
    Bnd *bd;
    std::list<Hdr*> h;
    uint32_t kct_l, bd_l;
    uint8_t kct_m, kcsndx_m, bd_m;
    std::map<uint32_t, Hdr*, Tid> hdr;
};

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
