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
#include "fa.h"

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

#define PENDING(k) ((*(k) & (ONE_CT - ONE_PENDING)) >> SMALL_SHFT)
#define ONE_PENDING (1ul << SMALL_SHFT)

#define STRAND_BIT (1ul << BIG_SHFT)
#define INDEX_MASK (STRAND_BIT - 1ul)

#define INFIOR_SHFT (BIG_SHFT + 1)
#define INFERIORITY (1ul << INFIOR_SHFT)
#define STRAND_POS (INFERIORITY - 1ul)
#define INFIOR_MASK (~STRAND_POS)
#define INFIOR INFERIORITY //((uint64_t)(iter+1) << INFIOR_SHFT)  //INFERIORITY
#define MAX_INFIOR 0xFFFFFF0000000000

#define B2POS_MASK 0x000000FFFFFFFFFF
#define KCT_B2POS(k) ((k)->fst & B2POS_MASK)

// XXX: Could use just one of these: not DISTINCT in non 1st iteration means MARKED.
#define MARKED 0x8000000000000000
#define DISTINCT 0x4000000000000000

#define ONE_CT (1ul << BIG_SHFT)
#define REMAIN(k) (((k)[1] >> BIG_SHFT) & 0x3FFFFF)

#define TSO_NT(k) ((k)[1] & B2POS_MASK)

#define IS_DISTINCT(k) (((k)[1] & DISTINCT) == DISTINCT)
#define ALL_SAME_NTS(k) (((k)[1] & MARKED) == MARKED)

#define _GET_NEXT_NT(kc, p) (((kc)->ts[((p) & B2POS_MASK)>>2] >> (((p) & 3) << 1)) & 3)


// Note: if PENDING(k) != 0 then we may or may not want to act.
#define IS_UQ(k) (REMAIN(k) == 1ul)

// assigns on purpose.
#define CONTINUED_UQ(k, kc) IS_UQ(k = kc->kct_scope[0])

#define FIRST_IN_SCOPE(k) (*(k) & (ONE_PENDING - 1ul))

// last in scope
#define LAST_IN_SCOPE(k) (FIRST_IN_SCOPE(k) + PENDING(k))

#define NT_OFFS(k) (TSO_NT(k) + LAST_IN_SCOPE(k))

#define IS_FIRST(k) ((*(k) & (ONE_PENDING - 1ul)) == 0ul)


#define IS_LAST(k) (REMAIN(k) == LAST_IN_SCOPE(k) && PENDING(k) == 0ul)


#define SAME_OR_UQ(k) (IS_UQ(k) || ALL_SAME_NTS(k))

// TODO:
#define PENDING_SAME(k) (IS_UQ(k) == false && ((k)[1] & (MARKED|DISTINCT)) == MARKED)

///////////////////////
//
#define DOWNSTREAM_ADJOINING(i, p) ((i) == (p))

#define FORMER_UQ_WAS_1ST(kc, i, p) (((*kc->bdit).e + i) == p)

#define IN_SCOPE_OF_LAST_BND(kc, p) (p < (*(kc)->bdit).s + (kc)->ext)

//FIXME: could pack dna in 32 bits for KEY_WIDTH 16.
packed_struct Mantra { // not yet covered by unique keys
    uint64_t dna; // dna after covered boundary
#ifdef DEBUG
    uint64_t end_dna; // dna at end of mantra (to check)
#endif
    uint32_t corr; // 'real' position correction
    uint32_t s, e; // start and end of mantra, position where we jump to.
};




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

#define read_b2(s, p) (((s)[(p)>>2] >> (((p) & 3) << 1)) & 3)

#define __rdsq(direction, b, s, b2pos, dna, rc) ({\
    b = b2pos;\
    ASSERT(b < s ## _l, return -EFAULT, ":%lu", b);\
    b = read_b2(s, b);\
    _seq_ ## direction (b, dna, rc);\
})

#define _get_new_kct(kc, k, dna, rc) ({\
    typeof(dna) __ndx, __t;\
    _get_ndx(__ndx, __t, dna, rc);\
    __t <<= BIG_SHFT - KEY_WIDTH; /*store orient and infior in t*/\
    k = kc->kct + kc->ndxkct[__ndx];\
    EPQ(dbgtsoffs == -2ul, "enter:[%lu]", TSO_NT(k));\
    ASSERT(REMAIN(k) != 0, res = -EFAULT; goto err);\
    ASSERT(IS_UQ(k) || LAST_IN_SCOPE(k) <= REMAIN(k), res = -EFAULT; goto err, "%lu, %lu", LAST_IN_SCOPE(k), REMAIN(k));\
    ASSERT(kc->ndxkct[__ndx] < kc->kct_l, res = -EFAULT; goto err);\
    EPQ0(TSO_NT(k) == dbgtsoffs, ">%s%s%s%s%sREMAIN:%lu\tNEXT_NT_NR:%lu(+1)\tPENDING:%lu\t",\
            IS_FIRST(k) ? "FIRST\t" : "", ALL_SAME_NTS(k) ? "ALL_SAME\t" : "", \
        IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",\
        IS_LAST(k) ? "LAST\t" : "",REMAIN(k), FIRST_IN_SCOPE(k), PENDING(k));\
    print_2dna(dna, rc, TSO_NT(k) == dbgtsoffs);\
    *k ^= (*k ^ __t) & STRAND_BIT;\
})

#define __rdndx(direction, ndx, b, s, b2pos, dna, rc) ({\
    dna = __rdsq(direction, b, s, b2pos, dna, rc);\
    _get_ndx(ndx, b, dna, rc);\
})

#define MARK_LAST(k)\
    if (IS_LAST(k)) {\
        if (IS_DISTINCT(k) == false)\
            k[1] |= MARKED;\
        EPQ(TSO_NT(k) == dbgtsoffs, " ==> dbgtsoffs was %sdistinct in last; pos reset",\
		k[1] & MARKED ? "_NOT_": "");\
        *k &= ~B2POS_MASK; /* reset position tracking */\
    }

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

packed_struct kct_ext {
    uint64_t m: 8;
    uint64_t l: 40;
    uint64_t* b2;
};

enum ensembl_parts {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR = 1};

struct Hdr {
    uint64_t s_s;
    uint32_t *part; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<Mantra> bnd; //
    uint32_t end_pos, end_corr, mapable; // mapable is only used in fa.cpp
    uint8_t p_l;
};

struct kct_t {
    char* id;
    uint8_t* ts; // next nts(nn)
    uint8_t* s; // all needed 2bit sequences in order (excluding Ns or first ones).
    uint32_t* ndxkct; // sparse array, complement independent index (ndx) => kct
    uint64_t* kct; // each 2 u64s with different usage in various stages, see below.
    uint64_t** kct_scope;
    uint64_t ts_l, s_l;
    uint32_t id_l, kct_l, uqct, pending;
    unsigned ext, iter; // not stored
    uint8_t id_m, s_m, ndxkct_m;
    uint8_t kct_m;
    std::list<Hdr*> h;
    std::list<Mantra>::iterator bdit;
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

int map_fq_se(struct seqb2_t*, char C*C);

// mantra.cpp
void show_mantras(kct_t C*C kc, Hdr *C h);
int insert_mantra(kct_t *C kc, Hdr* h);
void pot_mantra_end(kct_t *C kc, Hdr *C h, C uint64_t dna, C uint32_t b2pos);
#endif // RK_FA_H


