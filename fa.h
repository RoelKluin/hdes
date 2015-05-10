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
#include <errno.h> // ENOMEM
#include <list>
#include "seq.h"
#include "klib/khash.h"
#include "gz.h"

// length of the next NTs, when in sequence format, i.e. the lower bits
// contain a twobit rather than a index to a extended keycount (kct_ext)
#define B2LEN_OFFS_SHFT 59
#define B2LEN_OFFS (1ul << B2LEN_OFFS_SHFT)

// indicates whether (set) or not the kct contains a sequence or an index
#define FLAG_B2CT  (1ul << 58)

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
#define STRAND_SHFT 40
#define INDEX_MASK ((1ul << STRAND_SHFT) - 1ul)

// While some movement of next-NTs per key takes place - upwards
// movement of next-NTs that lie within range of unique indices and can
// therefore be skipped in future iterations, the number of next-Nts per
// key remain constant.


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

#define _get_ndx(ndx, wx, dna, rc) ({\
    wx = dna & KEYNT_STRAND;/* Store strand orientation. Central bit determines*/\
    ndx = wx ? dna : rc;    /* strand. Excise it out since its always the same */\
    ndx = ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
    ASSERT(ndx < KEYNT_BUFSZ, return -EFAULT, "0x%lx", ndx);\
    dbg = (ndx & INDEX_MASK) == dbgndx ? dbg | 4 : dbg & ~4;\
});

#define _update_kctndx(kc, ndx, wx, wi, dna, rc) ({\
    _get_ndx(ndx, wx, dna, rc);\
    wx <<= (STRAND_SHFT - KEY_WIDTH); /* strand orient. in position for storage*/\
    wx ^= (kc)->kctndx[ndx] & (1ul << STRAND_SHFT);\
    (kc)->kctndx[ndx] ^= wx; /* flip that bit */\
    wx = (kc)->kctndx[ndx];\
    wi = wx >> (STRAND_SHFT + 1); /* get this keys' inferiority */\
    ASSERT((wx & INDEX_MASK) < kc->kct_l, return -EFAULT, \
        "wx:%lx >= kc->kct_l:%lx, ndx:%lx, kctndx_m:%u, kctndx[ndx]:%lx, dna:%lx",\
        wx, kc->kct_l, ndx, (kc)->kctndx_m, (kc)->kctndx[ndx], dna);\
    wx &= INDEX_MASK; /* all else is index */\
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

#define ULL(x) ((unsigned __int128)(x))

#undef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#undef min
#define min(a,b) ((a) <= (b) ? (a) : (b))

#define DEBUG 1

packed_struct Bnd {
#ifdef DEBUG
    uint64_t at_dna; // dna before boundary (to check)
#endif
    uint64_t dna; // dna after boundary
    uint32_t s; // start
    int32_t l; // correction for position
    uint32_t i; // inferiority
};

packed_struct kct_ext {
    uint64_t m: 8;
    uint64_t l: 40;
    uint64_t* b2;
};

packed_struct Walker {
    uint32_t count;
    uint32_t excise_ct;
    uint32_t tmp_count;
};

enum ensembl_parts {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};

struct Hdr {
    uint8_t* s;
    uint32_t end_pos, s_l;
    uint32_t *part; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<uint32_t> bnd; //
    uint8_t p_l, s_m;
};

struct kct_t {
    Bnd* bd;
    char* id;
    uint8_t* ts; //later req
    uint64_t* kct; // different meaning later.
    kct_ext* kce; // early req
    Walker* wlkr; // later req
    uint64_t* wbuf; // later req
    uint64_t *kctndx;
    uint32_t bd_l, id_l, ext; // ext not stored
    uint32_t kce_l; //early req
    uint64_t ts_l, kct_l; // late req, continued req
    uint8_t kct_m, kce_m, kctndx_m, bd_m, id_m; // only bd_m is required, but not stored either
    std::list<Hdr*> h;
    std::list<uint32_t>::iterator bdit;
    // could be possible to move bnd here.
};
int fa_kc(kct_t*, void*, int (*) (void*), int (*) (int, void*));
int fa_index(seqb2_t*, uint32_t);
int write1(struct gzfh_t*, kct_t*);
int restore1(struct gzfh_t*, kct_t*);
int kct_convert(kct_t*);
#endif // RK_FA_H


