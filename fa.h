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
#include <queue>
#include "seq.h"
#include "klib/khash.h"
#include "gz.h"
#include "b6.h"

// length of the next NTs, when in sequence format, i.e. the lower bits
// contain a twobit rather than a index to a extended keycount (kct_ext)
//#define B2LEN_SHFT 58
//#define ONE_B2SEQ (1ul << B2LEN_SHFT)

//#define B2SEQ_MASK  (ONE_B2SEQ - 1ul)
//#define B2LEN_MASK  (~B2SEQ_MASK)

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
//#define BIG_SHFT 36
#define ORIENT_SHFT 36
#define SAME_CT_SHFT 37

// TODO: do not count, rather write pos on 1st occurance (1-based), in case of multiple, set length

#define DUP_BIT    0x8000000000000000 // Set if key is uniq for remaining sections (and pos stored).
#define STRAND_BIT 0x0000001000000000 // the bit to store the original orientation of ref once uniq.
#define MAX_INFIOR 0x7FFFFFD000000000 // inferiority per key - stays when unique
#define B2POS_MASK 0x0000000FFFFFFFFF // position, once unique

#if ((1 << ORIENT_SHFT) != STRAND_BIT)
#error "(1 << ORIENT_SHFT) != STRAND_BIT)"
#endif

// stored position is one-based to ensure a bit is set
#define NO_KCT -2u

#define B2POS_OF(k) ((k) & B2POS_MASK)

#define _PRNT_SEQ_BY_POS(kc, _p) ({\
    keyseq_t __seq = {0};\
    EPR0(Pfmt "\t", (pos_t)_p);\
    __seq.p = (pos_t)_p - KEY_WIDTH;\
    _build_key(kc, __seq, __seq.p, (pos_t)_p, __seq.t);\
    print_seq(&__seq);\
    -1;\
})

#define _PRNT_SEQ_BY_K(kc, k) ({\
    EPR0("koffs:%lu\t",  K_OFFS(kc, k));\
    _PRNT_SEQ_BY_POS(kc, B2POS_OF(*k));\
})

#define K_OFFS(kc, k) ((k) ? (k) - (kc)->kct : ~0ul)

// XXX ugly since returns with value, but for debugging..
#define _B2POS_OF(kc, k) ({\
        ASSERT(k - kc->kct < kc->kct_l, return -EFAULT);\
        ASSERT(k - kc->kct >= 0, return -EFAULT);\
        ASSERT(((*k) & B2POS_MASK) != NO_KCT, return -EFAULT);\
        if (B2POS_OF(*k) == dbgpos) {\
            EPR("observed dbggpos at koffs:%lu (%s:%u)", K_OFFS(kc, k), __FILE__, __LINE__);\
            _PRNT_SEQ_BY_POS(kc, dbgpos);\
        }\
        B2POS_OF(*k);\
})
//#define B2POS_lt(a, b) (_B2POS_OF(a) < _B2POS_OF(b))

// TODO: if a key is not unique, store pos or same seq and length in lower bits
//#define UQ_MASK    0x000000000000000F // How many same Nts occur
//#define SAME_SEQ   0X00000003FFFFFFF0

#define IS_DUP(k)      (*(k) & DUP_BIT)
#define IS_UQ(k)       (IS_DUP(k) == 0ul)
//#define IS_FIRST(k)    ((*(k) & B2POS_MASK) == 0ul)

#define INFIOR_SHFT (ORIENT_SHFT + 1)
#define INFERIORITY (1ul << INFIOR_SHFT)
#define STRAND_POS (INFERIORITY - 1ul)
#define INFIOR_MASK (~STRAND_POS)
#define INFIOR INFERIORITY

#define IS_DBG_K(kc, k) (*k == dbgk || K_OFFS(kc, k) == dbgndxkct)

#define DESCRIBE_KEY(kc, k, c) \
  EPR("[k:%lx, dbgk: %lx] %c\t%s", *k, K_OFFS(kc, k), c, IS_UQ(k) ? "UQ\t" : "")

// XXX: Could use just one of these: not DISTINCT in non 1st iteration means MARKED.
//#define MARKED 0x8000000000000000

//#define ALL_SAME_NTS(k) (IS_FIRST(k) && IS_DISTINCT(k))

//#define _GET_NEXT_NT(kc, p) (((kc)->ts[_B2POS_OF(p)>>2] >> (((p) & 3) << 1)) & 3)

//#define SAME_OR_UQ(k) (IS_UQ(k) || ALL_SAME_NTS(k))


packed_struct Mantra { // not yet covered by unique keys
    pos_t s, e; // start and end of mantra, position where we jump to.
    pos_t corr; // 'real' position correction
};

// While some movement of next-NTs per key takes place - upwards movement
// of next-NTs within range of unique indices and can therefore be skipped
// in future iterations - the number of next-Nts per key remains constant.

#define _seq_next(c, seq) ({\
    seq_t __b = c;\
    seq.rc = ((seq.rc << 2) & KEYNT_MASK) | (__b ^ 2);\
    (__b << KEYNT_TOP) | (seq.dna >> 2);\
})

#define _get_kct0(kc, seq, t, ndx, _act) ({\
    ndx = _get_ndx(t, seq.dna, seq.rc);\
    ASSERT(ndx < KEYNT_BUFSZ, print_seq(&seq); _act, Sfmt ", 0x%lx", ndx, KEYNT_BUFSZ);\
    dbg = (ndx == dbgndx || kc->ndxkct[ndx] == dbgndxkct ||\
            (kc->ndxkct[ndx] != NO_KCT && kc->kct[kc->ndxkct[ndx]] == dbgk)) ? dbg | 8 : dbg & ~8;\
    EPQ(dbg & 8, "observed dbg ndx " Sfmt " / ndxkct " Sfmt " / k 0x%lx: %s +%u",\
            ndx, kc->ndxkct[ndx], kc->kct[kc->ndxkct[ndx]], __FILE__, __LINE__);\
    ndx;\
})

#define _get_kct(kc, seq, t, _act) ({\
    seq_t __ndx;\
    __ndx = _get_kct0(kc, seq, t, __ndx, _act);\
    ASSERT(kc->ndxkct[__ndx] < kc->kct_l, print_seq(&seq); _act, Sfmt "\t" Sfmt, __ndx, kc->ndxkct[__ndx]);\
    kc->ndxkct + __ndx;\
})

#define ASSERT_SCOPE_RNG(kc, p, pend, _act)\
    do {\
        ASSERT(pend < (kc->s_l << 2), _act, Pfmt "/%lu?", pend, kc->s_l);\
        ASSERT(p < pend, _act, Pfmt " >= " Pfmt "?", p, pend);\
        ASSERT(p + kc->readlength - KEY_WIDTH >= pend, _act,\
                Pfmt " + %u < " Pfmt "?", p, kc->readlength - KEY_WIDTH, pend);\
    } while(0)

#define _build_key(kc, seq, p, pend, t)\
    ASSERT(p < pend, return -EFAULT, Pfmt " >= " Pfmt "?", p, pend);\
    ASSERT(pend < (kc->s_l << 2), return -EFAULT, Pfmt "/%lu?", pend, kc->s_l);\
    ASSERT(p >= 0, return -EFAULT, Pfmt " < 0 (" Pfmt ")?", p, pend);\
    do {\
        t = (kc->s[p>>2] >> ((p&3) << 1)) & 3;\
        seq.dna = _seq_next(t, seq);\
    } while (++p != pend)

#define _addtoseq(buf, b)\
    do {\
        if ((buf ## _l & 3) == 0) {\
            _buf_grow_err(buf, 1ul, 2, return -ENOMEM);\
            buf[buf ## _l>>2] = '\0';\
        }\
        buf[buf ## _l>>2] |= b << ((buf ## _l & 3) << 1);\
        ++(buf ## _l);\
    } while(0)

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
    uint32_t end_pos, end_corr; // mapable is only used in fa.cpp
    std::list<Mantra> *bnd; //
    uint8_t p_l;
};

struct HK {
    uint32_t kct;
    uint32_t hoffs;
    uint32_t ext;
};

struct kct_t {
    char* id;
    uint8_t* s; // all needed 2bit sequences in order (excluding Ns or first ones).
    seq_t* ndxkct; // somewhat sparse array, complement independent index (ndx) => kct
    uint64_t* kct; // each 2 u64s with different usage in various stages, see below.
    uint64_t** kct_scope;
    uint64_t s_l, totNts, h_l;
    uint32_t id_l, kct_l, hk_l, uqct, reeval, ext, last_uqct;
    unsigned readlength, iter;
    uint8_t id_m, s_m, ndxkct_m, h_m, kct_m, hk_m;
    Hdr* h;
    HK* hk;
    std::list<Mantra>::iterator bdit;
    // could be possible to move bnd here.
};

/* == kct in key_init stage: ==
 * key count in lowest bits.
 *
 * == unique boundary extension: ==
 * inferiority is needed, 23 bits originally. Highest bits for sorting?
 * one bit to mark unique
 * if unique, low bits contain position
 *
 * if not unique
 * key counts are needed.
 * required 2 bits to keep last nt, 1 bit to mark all same.
 *
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
int fa_read(struct gzfh_t*, kct_t*);
int fa_index(struct gzfh_t*, unsigned readlength);

int save_seqb2(struct gzfh_t*, kct_t*);
int load_seqb2(struct gzfh_t*, kct_t*);
int save_boundaries(struct gzfh_t*, kct_t*);
int load_boundaries(struct gzfh_t*, kct_t*);

int save_kc(struct gzfh_t*, kct_t*);
int load_kc(struct gzfh_t*, kct_t*);
int ammend_kc(struct gzfh_t*, kct_t*);

int map_fq_se(struct seqb2_t*, char C*C);

// mantra.cpp
int show_mantras(kct_t C*C kc, Hdr *C h);
int insert_mantra(kct_t *C kc, Hdr* h);
void pot_mantra_end(kct_t *C kc, Hdr *C h, C seq_t dna, C uint32_t b2pos);
#endif // RK_FA_H


