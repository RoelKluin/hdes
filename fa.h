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
#define STRAND_BIT (1ul << STRAND_SHFT)
#define INDEX_MASK (STRAND_BIT - 1ul)
#define INFIOR_SHFT (STRAND_SHFT + 1)
#define INFERIORITY (1ul << INFIOR_SHFT)
#define L_INFIOR (INFERIORITY - 1ul)
#define R_INFIOR (~(INFERIORITY - 1ul))

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

#define _incr_at_ndx(kc, left, ndx)\
    ASSERT(kc->wbuf[left-1] == ~0ul, return -EFAULT, "[%u/%u]", left-1, kc->ext);\
    kc->wbuf[left-1] = ndx;\
    ++kc->wlkr[ndx].tmp_count;

#define _extd_uq(kc, left, ndx)\
    left = kc->ext;\
    _incr_at_ndx(kc, left, ndx);

#define _get_ndx_and_strand(ndx, b, dna, rc) ({\
    b = dna & KEYNT_STRAND;/* Store strand orientation. Central bit determines*/\
    ndx = b ? dna : rc;    /* strand. Excise it out since its always the same */\
    ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
})

#define _get_ndx(ndx, wx, dna, rc) ({\
    ndx = _get_ndx_and_strand(ndx, wx, dna, rc);\
    ASSERT(ndx < KEYNT_BUFSZ, return -EFAULT, "0x%lx", ndx);\
    dbg = ((ndx == dbgndx) || ((kc)->kctndx[ndx] & INDEX_MASK) == dbgkctndx) ? dbg | 8 : dbg & ~8;\
    EPQ(dbg & 8, "observed dbgndx 0x%lx / dbgkctndx %lu", ndx, (kc)->kctndx[ndx] & INDEX_MASK);\
});

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

#define _verify_seq(b2pos, offs, ks, msg, action, ...) \
    if (kc->s) {\
        uint64_t p = b2pos + offs;\
        uint8_t sb2 = (kc->s[p>>2] >> ((p & 3) << 1)) & 3;\
        if (b2 != sb2) {\
            WARN("assertion 'b2 != sb2' failed " msg, __VA_ARGS__);\
            action;\
        }\
    }

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
    uint64_t dna; // dna after boundary, where we `jump' to
    uint32_t s;   // position at which we jump.
    uint32_t l;   // length of jump or real position correction.
    uint32_t corr;
};

packed_struct running {
    uint64_t infior;
    int left, uqct;
    unsigned ct;
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
    uint64_t s_s;
    uint32_t *part; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<uint32_t> bnd; //
#ifdef DEBUG
    uint32_t end_pos;
#endif
    uint8_t p_l;
};

struct kct_t {
    Bnd* bd;
    char* id;
    uint8_t* ts, *s; // next nts(nn), seqb2
    uint64_t* kct; // different meaning later.
    kct_ext* kce; // early req
    Walker* wlkr; // later req
    uint64_t* wbuf; // later req
    uint64_t *kctndx;
    uint32_t id_l, ext; // ext not stored
    uint32_t kce_l; //early req
    uint64_t bd_l, ts_l, s_l, kct_l; // late req, continued req
    uint8_t kct_m, kce_m, kctndx_m, bd_m, id_m, s_m; // only bd_m is required, but not stored either
    std::list<Hdr*> h;
    std::list<uint32_t>::iterator bdit;
    // could be possible to move bnd here.
};
void show_list(kct_t*, std::list<uint32_t> &bnd);
void free_kc(kct_t* kc);
int fa_read(struct seqb2_t*, kct_t*);
int fa_index(seqb2_t*);

int save_seqb2(struct gzfh_t*, kct_t*);
int load_seqb2(struct gzfh_t*, kct_t*);
int save_nextnts(struct gzfh_t*, kct_t*);
int load_nextnts(struct gzfh_t*, kct_t*);
int save_boundaries(struct gzfh_t*, kct_t*);
int load_boundaries(struct gzfh_t*, kct_t*);

int save_kc(struct gzfh_t*, kct_t*);
int load_kc(struct gzfh_t*, kct_t*);

int map_fq_se(struct seqb2_t*);

#endif // RK_FA_H


