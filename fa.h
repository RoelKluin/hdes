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
#include <string.h>
#include "seq.h"
//#include "klib/khash.h"
#include "gz.h"
#include "b6.h"

// TODO: INFIOR no longer needed and pos per contig - should fit in 32 bits

// in kct_convert(), all next NTs are put in a single buffer and kcts are
// converted to 2bit indices to their respective next Nts: Below this bit.
//#define BIG_SHFT 36

// stored position is one-based to ensure a bit is set
// this is the kct offset: can be 0.
#define NO_K -1u

#define K_OFFS(kc, k) ((k) - (kc)->kct)
//#define K_OFFS(kc, k) ((*(k) >> 32) & 0xffffffff)

#define IS_DUP(k) ({\
    decltype(*(k)) __t = *(k);\
    NB(_b2pos_of(__t) >= (KEY_WIDTH - 1) && __t != NO_K);\
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

#define get_kct(kc, seq, with_orient) ({\
    get_ndx(seq, with_orient);\
    NB(seq.t < KEYNT_BUFSZ);\
    kc->contxt_idx + seq.t;\
})

#define update_append_hkoffs(kc, e, hkoffs) ({\
    *hkoffs++ = K_OFFS(kc, e->tgtk);\
    unsigned __t = hkoffs - kc->hkoffs;\
    buf_grow_add(kc->hkoffs, 1ul, kc->kct_l);\
    kc->hkoffs + __t;\
})

struct Mantra { // not yet covered by unique keys
    uint32_t ho;   // which contig (header offfset)
    uint32_t s;    // start pos
    uint32_t corr; // 'real' position correction
    uint32_t e;    // end of mantra
};

// ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
enum ensembl_part {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR = 1};

packed_struct Hdr {
    uint32_t len, end; // 2bit length of this contig
    uint32_t ido;      // id contains offset to kc->id character for the header ID.
    uint32_t p_l:8;    // :4 How many parts in the ensembl format occurred (if one there's only an ID, format is unkown)
    uint32_t corr:24;  // total of N's + offsets, sums with with .len to contig lenght.
};

// variables only used while extending keys
struct Iter_t {
    uint32_t *tgtk; // pointer to where keys are stored that have not yet become unique.
                    // unique keys are initially removed and appended to kct, but after compression
                    // this is undone. There may be another strategy, swapping keys, but it's hard
                    // to implement because the contig needs to remain known.
    Mantra* obnd;   // for a next iteration old regions with mantras (kc->bnd) are stored here.
    uint8_t* s;     // offset of sequence for contig. first starts at 1st 2bit (not shifted).
    unsigned fk;    // first key of contig, only needed to recognise a same key within scope.
    unsigned obnd_l, obnd_m;
    unsigned moved; // offset to 1st of keys to be moved to start (the remaining non-uniques)
                    // this occurs unless keys are not in scope.
};

struct Key_t {
    Hdr* h;        // headers
    char* id;      // characters of headers

    uint8_t* s;    // all needed 2bit sequences in order (excluding Ns).
    uint32_t* contxt_idx; // strand independent k-mer to kct index.
                   // Later we may want to point to a combination in the non-observed for edits to
                   // this ndx or surroundings that do result in an ndx that does occur.
                   // non occurant are initally set to NO_K.
                   //
                   // If a key is found to be incorrect, One or more mismatches must have occurred
                   // within this key, + extension n. a 0-th (no) extension exists. If beyond
                   // readlength we cannot be conclusive about which position is right. (TODO:
                   // hash storage of remaining combinations of positions?)

    uint32_t* kct; // contains first occurrent position, strand and dupbit if not yet uniq.
                   // we could also store ndx for faster conversion, if 64 bit.

    uint32_t* ext_iter; // iterations per extension.
    uint32_t* hkoffs;   // kc->kct keys are kept ordered per contig. hkoffs indicates how many k's
                        // per extension per contig.

    Mantra* bnd;
    uint32_t h_l, id_l;
    uint64_t s_l, totNts;

    uint32_t kct_l, ext_iter_l, bnd_l, hkoffs_l;
    uint32_t ct, readlength;
    uint32_t ext, reserved;        // current to final extension.
    uint8_t contxt_idx_m, kct_m, ext_iter_m, bnd_m, hkoffs_m, h_m, id_m, s_m;
    // could be possible to move bnd here.
};

#define b2pos_of(k) ({\
    NB(k >= NT_WIDTH, "b2pos_of() %s:%i", __FILE__, __LINE__);\
    ((k) & B2POS_MASK);\
})

#define hdr_end_k(kc, h) ({\
    uint32_t __koffs = (kc)->hkoffs[(h) - (kc)->h];\
    NB(__koffs < (kc)->kct_l);\
    (kc)->kct + __koffs;\
})


#define build_ndx_kct(kc, seq, s, ...) ({\
    NB(seq.p >= NT_WIDTH && (b2pos_of(seq.p) >> 3) <= kc->s_l,\
            "build_ndx_kct() %s:%i\tb2pos:%x?", __FILE__, __LINE__, seq.p);\
    _build_ndx_kct(seq, s, ##__VA_ARGS__);\
})

void free_kc(Key_t* kc);
int fa_read(struct gzfh_t*, Key_t*);
int fa_index(struct gzfh_t*, uint64_t optm, unsigned readlength);

int save_seqb2(struct gzfh_t*, Key_t*);
int load_seqb2(struct gzfh_t*, Key_t*);

int save_kc(struct gzfh_t*, Key_t*);
int load_kc(struct gzfh_t*, Key_t*);

int map_fq_se(struct seqb2_t*, char C*C);

// mantra.cpp
int show_mantras(Key_t C*C kc, Mantra* obnd, unsigned obnd_l, Mantra* at);
void print_kct(Key_t *kc, Mantra* at, Iter_t* e, uint32_t* tk);
#endif // RK_FA_H


