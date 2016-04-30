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

// TODO: INFIOR no longer needed and pos per contig - should fit in 32 bits

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
//#define BIG_SHFT 36
#define ORIENT_SHFT 36
#define SAME_CT_SHFT 37

// TODO: do not count, rather write pos on 1st occurance (1-based), in case of multiple, set length

#define DUP_BIT    0x80000000 // Set if key is uniq for remaining sections (and pos stored).
#define STRAND_BIT 0x00000001 // the bit to store the original orientation of ref once uniq.
#define B2POS_MASK 0x7FFFFFFF // position, once unique

// stored position is one-based to ensure a bit is set
#define NO_KCT -2u

#define K_OFFS(kc, k) ((k) ? (k) - (kc)->kct : ~0ul)

#define IS_DUP(k) ({\
    decltype(*(k)) __t = *(k);\
    NB(b2pos_of(__t) > KEY_WIDTH - 1 && __t != NO_KCT);\
    __t & DUP_BIT;\
})
#define IS_UQ(k)       (IS_DUP(k) == 0ul)

#define IS_DBG_K(kc, k) (*k == dbgk || K_OFFS(kc, k) == dbgndxkct)

#define DESCRIBE_KEY(kc, k, c) \
  EPR("[k:%lx, dbgk: %lx] %c\t%s", *k, K_OFFS(kc, k), c, IS_UQ(k) ? "UQ\t" : "")

packed_struct Mantra { // not yet covered by unique keys
    pos_t s, e; // start and end of mantra, position where we jump to.
    pos_t corr; // 'real' position correction
};

// While some movement of next-NTs per key takes place - upwards movement
// of next-NTs within range of unique indices and can therefore be skipped
// in future iterations - the number of next-Nts per key remains constant.

static inline void
seq_next(struct keyseq_t &seq)
{
    seq.rc = ((seq.rc << 2) & KEYNT_MASK) | (seq.t ^ 2);
    seq.dna = seq.t << KEYNT_TOP | seq.dna >> 2;
}

#define get_kct0(kc, seq, ndx) ({\
    ndx = get_ndx(seq);\
    NB(ndx < KEYNT_BUFSZ);\
    ndx;\
})

#define get_kct(kc, seq) ({\
    seq_t __ndx;\
    __ndx = get_kct0(kc, seq, __ndx);\
    NB(kc->ndxkct[__ndx] < kc->kct_l);\
    kc->ndxkct + __ndx;\
})

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

#define _prev_or_bnd_start(b) (b.prev ? b2pos_of(*b.prev) : (*b.it).s + KEY_WIDTH - 1)

#define in_scope(kc, fst, nxt) ((nxt) - (fst) - 1u < (kc)->ext)

enum ensembl_parts {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR = 1};

struct Hdr {
    uint64_t s_s;
    uint32_t *part; //ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
    std::list<Mantra> *bnd; // XXX: waarom van header? Mantra *end in plaats, list in kct_t?
    Mantra* mend;
    uint32_t end_pos; // mapable is only used in fa.cpp
    uint8_t p_l;
};

// to look up key offset to 1) respective header and 2) for what extension the keys became uniq
struct HK {
    uint32_t hoffs;
    uint32_t ext;
    uint32_t koffs;
};

struct Bnd {
    pos_t *sk;
    Hdr *h;
    pos_t *prev;
    std::list<Mantra>::iterator it;
};

struct kct_t {
    char* id;
    uint8_t* s; // all needed 2bit sequences in order (excluding Ns or first ones).
    seq_t* ndxkct; // somewhat sparse array, complement independent index (ndx) => kct
    pos_t* kct;
    uint64_t s_l, totNts;
    uint32_t id_l, kct_l, hk_l, h_l, uqct, reeval, ext, last_uqct;
    unsigned readlength, iter;
    uint8_t id_m, s_m, ndxkct_m, h_m, kct_m, hk_m;
    Hdr* h;
    HK* hk;
    std::list<Mantra> bnd;
    // could be possible to move bnd here.
};

// TODO: pos rshift may not be necessary.
static inline pos_t
b2pos_of(kct_t C*C kc, pos_t C*C k)
{
    NB(k - kc->kct < kc->kct_l);
    NB(k - kc->kct >= 0);
    NB(*k != NO_KCT);
    
    return (*k & B2POS_MASK) >> 1;
}

static inline pos_t
b2pos_of(pos_t C k)
{
    NB(k != NO_KCT);
    return (k & B2POS_MASK) >> 1;
}

static seq_t
build_ndx_kct(keyseq_t &seq, uint8_t const*const s)
{
    pos_t p = seq.p;
    seq.p -= KEY_WIDTH;
    build_key(s, seq, p);
    seq_t ndx = get_ndx(seq);
    NB(ndx < KEYNT_BUFSZ);
    return ndx;
}

void free_kc(kct_t* kc);
int fa_read(struct gzfh_t*, kct_t*);
int fa_index(struct gzfh_t*, uint64_t optm, unsigned readlength);

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


