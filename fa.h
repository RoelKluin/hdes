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
#define B2LEN_SHFT 58
#define ONE_B2SEQ (1ul << B2LEN_SHFT)

#define B2SEQ_MASK  (ONE_B2SEQ - 1ul)
#define B2LEN_MASK  (~B2SEQ_MASK)

// indicates whether (set) or not the kct contains a sequence or an index
#define FLAG_B2CT  (1ul << 58)

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
#define BIG_SHFT 40
#define SMALL_SHFT 24

#define STRAND_BIT (1ul << BIG_SHFT)
#define INDEX_MASK (STRAND_BIT - 1ul)

#define INFIOR_SHFT (BIG_SHFT + 1)
#define INFERIORITY (1ul << INFIOR_SHFT)
#define STRAND_POS (INFERIORITY - 1ul)
#define INFIOR_MASK (~STRAND_POS)
#define INFIOR INFERIORITY //((uint64_t)(iter+1) << INFIOR_SHFT)  //INFERIORITY
#define MAX_INFIOR 0xFFFFFF0000000000

#define ONE_CT (1ul << BIG_SHFT)
#define B2POS_MASK 0x000000FFFFFFFFFF
#define KCT_B2POS(k) ((k)->fst & B2POS_MASK)

#define MARKED 0x8000000000000000
#define DISTINCT 0x4000000000000000

#define TSO_NT(k) ((k)[1] & B2POS_MASK)

#define IS_DISTINCT(k) (((k)[1] & DISTINCT) == DISTINCT)
#define ALL_SAME_NTS(k) (((k)[1] & MARKED) == MARKED)

#define _GET_NEXT_NT(kc, p) (((kc)->ts[((p) & B2POS_MASK)>>2] >> (((p)&3) << 1)) & 3)

//WAS_UQ is not necessary. Decrement no longer occurs if unique
#define IS_UQ(k) (((k)[1] & ~(MARKED|DISTINCT|ONE_CT)) <= ONE_CT)

// N.B. one baased or zero based dependent on whther it is used before or after increment/decrement
#define REMAIN(k) (((k)[1] & (DISTINCT-ONE_CT)) >> BIG_SHFT)

#define NEXT_NT_NR(k) (*(k) & B2POS_MASK)

#define IS_FIRST(k) (NEXT_NT_NR(k) == 0ul)

#define IS_LAST(k) (REMAIN(k) == NEXT_NT_NR(k))
#define WAS_LAST(k) (REMAIN(k) == NEXT_NT_NR(k) +1)


#define SAME_OR_UQ(k) (IS_UQ(k) || ALL_SAME_NTS(k))

// TODO:
#define PENDING_SAME(k) (IS_UQ(k) == false && ((k)[1] & (MARKED|DISTINCT)) == MARKED)

// While some movement of next-NTs per key takes place - upwards movement
// of next-NTs within range of unique indices and can therefore be skipped
// in future iterations - the number of next-Nts per key remains constant.


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

// if kc->ext, rotate to zero
#define KC_ROT_NEXT(kc, x) ((x+1) & -((x+1) != (kc)->ext))
#define KC_ROT(kc, x) (x &= -(++x != (kc)->ext))
#define KC_LEFT_ROT_NEXT(kc, x) (x + (-(x == 0) & (kc)->ext) - 1)
#define KC_LEFT_ROT(kc, x) (x = KC_LEFT_ROT_NEXT(kc, x))


#define _get_new_kct(kc, k, dna, rc) ({\
    typeof(dna) __ndx, __t;\
    _get_ndx(__ndx, __t, dna, rc);\
    __t <<= BIG_SHFT - KEY_WIDTH; /*store orient and infior in t*/\
    k = kc->kct + kc->ndxkct[__ndx];\
    ASSERT(REMAIN(k) != 0, return -EFAULT, "%lu", k[1] & B2POS_MASK);\
    ASSERT(IS_UQ(k) || NEXT_NT_NR(k) <= REMAIN(k), return -EFAULT, "%lu", k[1] & B2POS_MASK);\
    ASSERT(kc->ndxkct[__ndx] < kc->kct_l, return -EFAULT);\
    EPQ0(TSO_NT(k) == dbgtsoffs, "%s%s%s%s%sREMAIN:%lu\tNEXT_NT_NR:%lu(+1)\t",\
            IS_FIRST(k) ? "FIRST\t" : "", ALL_SAME_NTS(k) ? "ALL_SAME\t" : "", \
        IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",\
        IS_LAST(k) ? "LAST\t" : "",REMAIN(k), NEXT_NT_NR(k));\
    print_2dna(dna, rc, TSO_NT(k) == dbgtsoffs);\
    *k ^= (*k ^ __t) & STRAND_BIT;\
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

//KHASH_MAP_INIT_INT64(UQCT, unsigned)

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
    uint32_t l;   // length of jump.
    uint32_t corr;// real position correction
};

packed_struct kct_ext {
    uint64_t m: 8;
    uint64_t l: 40;
    uint64_t* b2;
};

enum ensembl_parts {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};

struct Hdr {
    uint64_t s_s;
    uint32_t *part; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<uint32_t> bnd; //
    uint32_t end_pos, mapable; // mapable is only used in fa.cpp
    uint8_t p_l;
};

struct kct_t {
    Bnd* bd;
    char* id;
    uint8_t* ts; // next nts(nn)
    uint8_t* s; // all needed 2bit sequences in order (excluding Ns or first ones).
    uint32_t* ndxkct; // sparse array, complement independent index (ndx) => kct
    uint64_t* kct; // each 2 u64s with different usage in various stages, see below.
    uint64_t** kct_scope;
    uint64_t ts_l, s_l;
    uint32_t id_l, bd_l, kct_l, uqct, pending;
    unsigned ext; // not stored
    uint8_t bd_m, id_m, s_m, ndxkct_m; // only bd_m is required, but not stored either
    uint8_t kct_m;
    std::list<Hdr*> h;
    std::list<uint32_t>::iterator bdit;
    // could be possible to move bnd here.
};
/* == kct in key_init stage: ==
 * next_nt (sequence) is stored in first u64 and 2nd up to 6 highest bits.
 * 6 highest bits contain the count. if more then 61 Nts, all 6 highest are set,
 * count is moved to lowest bits of 2nd u64. The 1st u64 receives the address of a
 * pointer to u64 that receives the storage and grows dynamically in powers of 2.
 * its size (2-to the powers) is also stored in 2nd u64, 32th bit and above.
 *
 * == kct conversion: ==
 * During conversion all nextNts are concatenated in kc->ts. The 2nd u64 receives
 * the kc->ts offset in its low 40 bits, its high 24 bits receive the nextNt count.
 * The 1st u64 of the kct receives the complement independent index (ndx). This is
 * just temporary - it allows more efficient ndxkct storage.
 *
 * == unique boundary extension: ==
 * The 2nd u64 of kct contains kc->ts offset in its low 40, nextNt count in high 24.
 * In the scope of another unique key, remaining non unique positions are no longer
 * considered. They can be skipped and therefore their respective nextNts are
 * excised and moved upwards. The nextNt count is decremented.
 *
 * The 1st u64 is unset and will be used for keeping track which nextnts are passed
 * for this index in the low 40 bits. The inferiority is stored in the highest 23th
 * bit and the 40th bit contains the orientation of the complement independent index
 * on reference once unique. Once unique the low 40 bits are used instead for the
 * (genomic) b2pos storage.
 * [0]: pos:40, strand:1, infior:23;
 * [1] ts_offs:40, remain:24;
 */
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
int ammend_kc(struct gzfh_t*, kct_t*);

int map_fq_se(struct seqb2_t*);

#endif // RK_FA_H


