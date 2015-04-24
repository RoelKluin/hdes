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

#define __seq_le_n(b, dna, rc) ({\
    rc = ((b ^ 2) << KEYNT_TOP) | (rc >> 2);\
    ((dna << 2) & KEYNT_MASK) | b;\
})
#define __seq_le_p(b, dna, rc) ({\
    rc = ((rc << 2) & KEYNT_MASK) | (b ^ 2);\
    (b << KEYNT_TOP) | (dna >> 2);\
})

#ifdef B2_LITTLE_ENDIAN
# define _seq_next __seq_le_n
# define _seq_prev __seq_le_p
#else
# define _seq_next __seq_le_p
# define _seq_prev __seq_le_n
#endif


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

#define _get_ndx_and_strand(ndx, b, dna, rc) ({\
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
        /**/if (dbg > 1) { fputc(b6(b << 1), stderr); }/**/\
    }\
    /**/if (dbg > 1) { fputc('\n', stderr); }/**/\
})
#define _addtoseq(buf, b)\
    if ((buf ## _l & 3) == 0) {\
        _buf_grow_err(buf, 1ul, 2, return -ENOMEM);\
        buf[buf ## _l>>2] = '\0';\
    }\
    buf[buf ## _l>>2] |= b << ((buf ## _l & 3) << 1);\
    ++(buf ## _l)

KHASH_MAP_INIT_INT64(UQCT, unsigned)

#define UNINITIALIZED (~0u)

#define INFERIORITY_BIT (1u<<31)
#define ORIENTATION (1u<<30)

#define ULL(x) ((unsigned __int128)(x))

#define NR_MASK 0x3fffffff
#define TYPE_MASK 0xc0000000
#undef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#undef min
#define min(a,b) ((a) <= (b) ? (a) : (b))

#define DEBUG 1

#ifdef DEBUG
# define DEBUG_ASSIGN_ENDDNA(at_dna_ref, dna) (at_dna_ref) = (dna)
# define ASSIGN_BD(_bd, _t, _d, _s, _l, _i, _e) \
        (_bd) = {.t = _t, .dna = _d, .s = _s, .l = _l, .i = _i, .at_dna = _e};

# define _getxtdndx0(kc, ndx) (kc)->kcsndx[(ndx)]
# define get_w(wlkr, kc, ndx) (wlkr)[(kc)->kcsndx[(ndx)]]
# define _get_kct(kc, ndx) \
        /*ASSERT((kc)->kcsndx[(ndx)] < (kc)->kct_l, return -EINVAL);*/\
        (kc)->kct[(kc)->kcsndx[(ndx)]]
# define _getxtdndx(kc, ndx, dna, rc) ({\
        ndx = _get_ndx((ndx), (dna), (rc));\
        /*EPQ(ndx == 0x47ef9, "%s:%u <========\n", hdr, pos);*/\
        ASSERT(ndx < (1ul << KEYNT_BUFSZ_SHFT), return -EINVAL);\
        ndx;\
})
#else
# define DEBUG_ASSIGN_ENDDNA(at_dna_ref, dna) //nothing
#define ASSIGN_BD(_bd, _t, _d, _s, _l, _i, _e) \
        (_bd) = {.t = _t, .dna = _d, .s = _s, .l = _l, .i = _i};

# define _getxtdndx0(kc, ndx) (ndx)
# define get_w(wlkr, kc, ndx) (wlkr)[(ndx)]
# define _get_kct(kc, ndx) (kc)->kct[(ndx)]
# define _getxtdndx(kc, ndx, dna, rc) ({\
        ndx = _get_ndx(ndx, dna, rc);\
        ASSERT(ndx < (1ul << KEYNT_BUFSZ_SHFT), return -EINVAL);\
        (kc)->kcsndx[ndx];\
})

#endif

typedef packed_struct bnd_t {
    uint64_t t: 8; // type
    uint64_t dna: 56;
    uint32_t s; // start
    uint32_t l; // length
    uint32_t i; // inferiority
#ifdef DEBUG
    uint64_t at_dna;
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
};

packed_struct Walker {
    uint32_t count;
    uint32_t infior: 22; //b2pos start and end of range
    uint32_t tmp_count: 10;
};

struct Hdr {
    uint8_t* s;
    uint32_t end_pos, s_l;
    uint16_t part[10]; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<uint32_t> bnd; //
    std::list<uint32_t>::iterator bdit;
    uint8_t p_l, s_m;
};

struct Tid {
    bool operator () (char*a, char* b) const 
    {
        return strcmp(a, b) < 0;
    }
};

struct kct_t {
    Bnd* bd;
    char* id;
    Kct* kct;
    Walker* wlkr;
    uint32_t* wbuf;
    uint32_t *kcsndx;
    uint32_t bd_l, id_l, kct_l, ext;
    uint8_t kct_m, kcsndx_m, bd_m, id_m;
    Tid tid;
    std::list<Hdr*> h;
    std::map<char*, Hdr*, Tid> hdr;
};

int fa_index(seqb2_t *seq, uint32_t blocksize);
#endif // RK_FA_H


