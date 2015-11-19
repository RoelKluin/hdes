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

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
#define BIG_SHFT 34

// TODO: do not count, rather write pos on 1st occurance (1-based), in case of multiple, set length

#define DUP_BIT    0x8000000000000000 // Set if key is uniq for remaining sections (and pos stored).
#define STRAND_BIT 0x0000000400000000 // the bit to store the original orientation of ref once uniq.
#define MAX_INFIOR 0x7FFFFFF800000000 // inferiority per key - stays when unique
#define B2POS_MASK 0x00000003FFFFFFFF // position, once unique

// TODO: if a key is not unique, store pos or same seq and length in lower bits
#define UQ_MASK    0x000000000000000F // How many same Nts occur
#define SAME_SEQ   0X00000003FFFFFFF0

#define IS_DUP(k)      (*(k) & DUP_BIT)
#define IS_UQ(k)       (IS_DUP(k) == 0ul)
#define IS_FIRST(k)    ((*(k) & B2POS_MASK) == 0ul)

#define INFIOR_SHFT (BIG_SHFT + 1)
#define INFERIORITY (1ul << INFIOR_SHFT)
#define STRAND_POS (INFERIORITY - 1ul)
#define INFIOR_MASK (~STRAND_POS)
#define INFIOR INFERIORITY

#define K_OFFS(kc, k) ((k) ? (k) - (kc)->kct : ~0ul)
#define IS_DBG_K(kc, k) (K_OFFS(kc, k) == dbgk)

#define NO_KCT -2u

#define DESCRIBE_KEY(kc, k, c) \
  EPR("[k:%lu] %c\t%s", K_OFFS(kc, k), c, IS_FIRST(k) || IS_UQ(k) ? "UQ\t" : "")

// XXX: Could use just one of these: not DISTINCT in non 1st iteration means MARKED.
#define MARKED 0x8000000000000000

#define ALL_SAME_NTS(k) (IS_FIRST(k) && IS_DISTINCT(k))

#define _GET_NEXT_NT(kc, p) (((kc)->ts[((p) & B2POS_MASK)>>2] >> (((p) & 3) << 1)) & 3)

#define SAME_OR_UQ(k) (IS_UQ(k) || ALL_SAME_NTS(k))


//FIXME: could pack dna in 32 bits for KEY_WIDTH 16.
packed_struct Mantra { // not yet covered by unique keys
    uint32_t s, e; // start and end of mantra, position where we jump to.
    uint32_t corr; // 'real' position correction
};




// While some movement of next-NTs per key takes place - upwards movement
// of next-NTs within range of unique indices and can therefore be skipped
// in future iterations - the number of next-Nts per key remains constant.


#define __seq_le_n(c, dna, rc) ({\
    typeof(dna) __b = c;\
    rc = ((__b ^ 2) << KEYNT_TOP) | (rc >> 2);\
    ((dna << 2) & KEYNT_MASK) | __b;\
})
#define __seq_le_p(c, dna, rc) ({\
    typeof(dna) __b = c;\
    rc = ((rc << 2) & KEYNT_MASK) | (__b ^ 2);\
    (__b << KEYNT_TOP) | (dna >> 2);\
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

#define _get_new_kct(kc, k, t, dna, rc) ({\
    typeof(dna) __ndx, __t;\
    _get_ndx(__ndx, __t, dna, rc);\
    ASSERT(kc->ndxkct[__ndx] < kc->kct_l, res = -EFAULT; goto err);\
    k = kc->kct + kc->ndxkct[__ndx];\
    if (IS_DBG_K(kc, k)) {\
        DESCRIBE_KEY(kc, k, '>');\
        print_2dna(dna, rc, IS_DBG_K(kc, k));\
    }\
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
    uint32_t end_pos, end_corr, mapable, total; // mapable is only used in fa.cpp
    uint8_t p_l;
};

struct kct_t {
    char* id;
    uint8_t* s; // all needed 2bit sequences in order (excluding Ns or first ones).
    uint32_t* ndxkct; // sparse array, complement independent index (ndx) => kct
    uint64_t* kct; // each 2 u64s with different usage in various stages, see below.
    uint64_t** kct_scope;
    uint64_t s_l;
    uint32_t id_l, kct_l, uqct, pending;
    unsigned readlength, iter;
    uint8_t id_m, s_m, ndxkct_m;
    uint8_t kct_m;
    std::list<Hdr*> h;
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
void show_mantras(kct_t C*C kc, Hdr *C h);
int insert_mantra(kct_t *C kc, Hdr* h);
void pot_mantra_end(kct_t *C kc, Hdr *C h, C uint64_t dna, C uint32_t b2pos);
#endif // RK_FA_H


