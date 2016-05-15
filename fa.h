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
// this is the kct offset: can be 0.
#define NO_KCT -1u

#define K_OFFS(kc, k) ((k) ? (k) - (kc)->kct : ~0ul)

#define IS_DUP(k) ({\
    decltype(*(k)) __t = *(k);\
    NB(b2pos_of(__t) >= (KEY_WIDTH - 1) && __t != NO_KCT);\
    __t & DUP_BIT;\
})
#define IS_UQ(k)       (IS_DUP(k) == 0ul)

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
    uint32_t __ndx;\
    __ndx = get_kct0(kc, seq, __ndx);\
    NB(kc->contxt_idx[__ndx] < kc->kct_l);\
    kc->contxt_idx + __ndx;\
})

#define _addtoseq(buf, seq)\
    do {\
        seq_next(seq);\
        if ((buf ## _l & 3) == 0) {\
            _buf_grow_err(buf, 1ul, 2, return -ENOMEM);\
            buf[buf ## _l>>2] = '\0';\
        }\
        buf[buf ## _l>>2] |= seq.t << ((buf ## _l & 3) << 1);\
    } while(0)

#define next_seqpos(buf, seq)\
    do {\
        seq.p += 2;\
        ++(buf ## _l);\
    } while(0)

#define DEBUG 1

#define _prev_or_bnd_start(b) (b.prev ? b2pos_of(*b.prev) : (*b.it).s + KEY_WIDTH - 1)

#define get_kend(kc) (kc->kct + kc->kct_l - kc->last_uqct)

#define in_scope(kc, fst, nxt) ({\
    uint32_t __f = fst, __n = nxt;\
    NB(__f < __n);\
    NB(__n <= kc->s_l);\
    (nxt) - (fst) - 1u < (kc)->extension;\
})

packed_struct Mantra { // not yet covered by unique keys
    uint32_t s; // start pos
    uint32_t corr; // 'real' position correction
    uint32_t e; // end pos
    uint32_t ke; // offset to end k of mantra
};

/* TODO:
// a region of sequence, not yet covered by unique keys
packed_struct Mantra {
    uint64_t s, e;  // 2bit start and end positions of a mantra.
};

// section of sequence without ambiguous nucleotides
struct Assembly
{
    std::list<Mantra> mantra;
    uint64_t corr;  // correction needed to get from 2bit position to 'real' position.
};*/

// ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
enum ensembl_part {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR = 1};

struct Hdr {
    uint32_t ido;     // id contains offset to kc->id character for the header ID.
    uint8_t p_l;      // How many parts in the ensembl format occurred (if one there's only an ID, format is unkown)
};

/*
 * the kc->kct keys are ordered. Initially upon occurance on the genome
 */
// to look up key offset to 1) respective header and 2) for what extension the keys became uniq
struct HK {
    uint32_t hoffs;
    uint32_t ext;
    uint32_t len; // 2bit offset for this contig
    uint32_t koffs;
};

struct Bnd {
    uint32_t *sk;
    uint8_t* s;
    uint32_t *prev;
    int rot;
    std::list<Mantra>::iterator it;
};

struct kct_t {
    char* id;      // characters of headers
    uint8_t* s;    // all needed 2bit sequences in order (excluding Ns or first ones).
    uint32_t* contxt_idx;
                   // Later we may want to point to a combination in the non-observed for edits to
                   // this ndx or surroundings that does result in an ndx that does occur.

    uint32_t* kct;
                   // non occurant are initally set to NO_KCT. see Extension below.

    uint64_t* ext; // ndx offset for key extensions. If a position is after this u64 of the nth
                   // extension n u64 (but before the next, or NO_KCT) then the key became unique
                   // with this extension. If a key is found to be incorrect, One or more
                   // mismatches must have occurred within this key, + extension n.
                   // a 0-th (no) extension exists. If beyond readlength we cannot be conclusive.
    uint64_t s_l, totNts;
    uint32_t id_l, kct_l, hk_l, h_l, uqct, last_uqct;
    unsigned readlength, iter, ext_l, ext_m, extension;
    uint8_t id_m, s_m, contxt_idx_m, h_m, kct_m, hk_m;
    Hdr* h;
    HK* hk;
    std::list<Mantra>* bnd;
    // could be possible to move bnd here.
};

// TODO: pos rshift may not be necessary.
static inline uint32_t
b2pos_of(kct_t C*C kc, uint32_t C*C k)
{
    NB(k - kc->kct < kc->kct_l);
    NB(k - kc->kct >= 0);
    NB(*k != NO_KCT);

    return (*k & B2POS_MASK) >> 1;
}

static inline uint32_t
b2pos_of(uint32_t C k)
{
    NB(k != NO_KCT);
    return (k & B2POS_MASK) >> 1;
}

static uint32_t
build_ndx_kct(keyseq_t &seq, uint8_t const*const s)
{
    uint32_t p = seq.p;
    seq.p -= KEY_WIDTH;
    build_key(s, seq, p);
    uint32_t ndx = get_ndx(seq);
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
void pot_mantra_end(kct_t *C kc, Hdr *C h, C uint32_t dna, C uint32_t b2pos);
#endif // RK_FA_H


