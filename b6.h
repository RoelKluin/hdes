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
#ifndef RK_B6_H
# define RK_B6_H
#include <stdint.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

# define BITS_PER_LONG 64

# define B6_MASK			0x6

# define B6_UC				0x40
# define B6_LC				0x60

# define B6_DNA				0x10
# define B6_RNA				0x11

# define B6_A2B				0x75
# define B6_B2A				0x35

# define B6_ALT_CASE			0x20

// XXX TODO XXX: if we do not search for a max, but rather select a key based on
// keycount lookup, then the cutoff is not needed and can be removed. However...

// KEY_WIDTH must be odd - or 2nd bit of central Nt is not always flipped in its
// complement - the alternative is the twisted halfdev conversion, but this is cheaper
#ifndef KEY_LENGTH
# define KEY_LENGTH 15                      // <= this many Nts are used as key
#endif

#if KEY_LENGTH > 16
#error "KEY_LENGTH > 16 not supported"
#endif

#define DEVIANT_BIT_MASK 0xaaaaaaaaU
#define ODD_BYTE_MASK 0x00ff00ffU
#define ODD_NIBBLE_MASK 0x0f0f0f0fU
#define ODD_TWOBIT_MASK 0x33333333U

#define DUP_BIT    0x80000000 // Set if key is uniq for remaining sections (and pos stored).
#define STRAND_BIT 0x00000001 // the bit to store the original orientation of ref once uniq.
#define B2POS_MASK 0x7FFFFFFE // position, once unique

#define KEY_CUTOFF 0                        // cut off, to get random keys, [XXX: keep for assembly]
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF) // entire key maximized on (after conversion)
#define NT_WIDTH (KEY_WIDTH << 1)
#define KEYNT_TOP (NT_WIDTH - 2)

// maximum length of the read name
#define SEQ_MAX_NAME_ETC    (1u << 14)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ_SHFT ((KEY_LENGTH * 2) - 1)
#define KEYNT_BUFSZ (1ul << KEYNT_BUFSZ_SHFT)
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1ul) // also contxt_idx max.

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_MASK ((1ul << NT_WIDTH) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)
#define KEYNT_TRUNC_UPPER (~HALF_KEYNT_MASK & KEYNT_TRUNC_MASK)


#define SNDX_TRUNC_MASK (KEYNT_TRUNC_MASK | KEYNT_BUFSZ)
#define HALF_SNDX_MASK (HALF_KEYNT_MASK | KEYNT_STRAND)
#define SNDX_TRUNC_UPPER (~HALF_SNDX_MASK & SNDX_TRUNC_MASK)

#define _b2pos_of(pos) ((pos) & B2POS_MASK)


// set b to 2bit for c, returns true if twobit
#define B6(b, c) ({\
    b = c ^ ((c | B6_UC) & B6_LC); /* flip case bits */\
    b ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);\
    b ^= -((c & 0x8e) == 0x4) & B6_RNA; /* b is now converted to or from 2bit */\
    (b | B6_MASK) == B6_MASK;/*test whether the result is a twobit*/\
})


// zero [min] for N - but also for A.
//
#define b6N0(t, c) ({\
        t = c;\
        t ^= ((c | B6_UC) & B6_LC);/* FIXME maybe simplify with line below?*/\
        t ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);\
        t ^= -((c & 0x8e) == 0x4) & B6_RNA;\
        -((t | B6_MASK) == B6_MASK) & (t >> 1);\
})

#define add_b6N0(t, c, dna, rev) ({\
        t = c;\
        t ^= ((c | B6_UC) & B6_LC);\
        t ^= -((c & 0x8e) == 0x4) & B6_RNA;\
        t ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);\
        t = -((t | B6_MASK) == B6_MASK) & (t >> 1);\
        dna = dna << 2 | t;\
        rev = (uint32_t)t << KEYNT_TOP | rev >> 2;\
})

#define add_b(t, dna, rev) ({\
        dna = dna << 2 | t;\
        rev = (uint32_t)t << KEYNT_TOP | rev >> 2;\
})

# define ischar_ignore_cs(c, c2) (((c ^ c2) & ~B6_ALT_CASE) == 0)

#define REVSEQ(dna, m) ({\
    m = ODD_TWOBIT_MASK;\
    dna = (dna & m) << 2 | (dna & ~m) >> 2;\
    m = ODD_NIBBLE_MASK;\
    dna = (dna & m) << 4 | (dna & ~m) >> 4;\
    m = ODD_BYTE_MASK;\
    asm ("bswap %0" : "=r" (dna) : "0" (dna));\
    dna >> ((sizeof(dna) << 3) - NT_WIDTH);\
})

inline uint32_t revcmp(uint32_t dna) /* Reverse Complement, is ok. */
{
    uint32_t m;
    dna ^= DEVIANT_BIT_MASK;
    return REVSEQ(dna, m);
}

struct __attribute__ ((__packed__)) keyseq_t {
    uint32_t p;
    uint32_t dna, rc, t;
};

static inline void
get_next_nt_seq(uint8_t const*const s, struct keyseq_t &seq)
{
    if (seq.p & 1) {
        fprintf(stderr, "get_next_nt_seq(): orient already set\n");
        fflush(stderr);
        raise(SIGTRAP);
    }
    seq.t = (s[seq.p>>3] >> (seq.p&6)) & 3;
    seq.rc = ((seq.rc << 2) & KEYNT_MASK) | (seq.t ^ 2);
    seq.dna = seq.t << KEYNT_TOP | seq.dna >> 2;
}

extern inline void get_next_nt_seq(uint8_t const*const s, struct keyseq_t &seq);
unsigned b6(unsigned c);
unsigned b6_spec(unsigned c, unsigned cs, unsigned no_u);

#define build_seq(s, seq, pend)\
    do {\
        seq.p -= NT_WIDTH;\
        do {\
            get_next_nt_seq(s, seq);\
            seq.p += 2;\
        } while (seq.p != pend);\
    } while (0)


// if with_orient is 1, the orientation is stored in the first bit of seq.p (
#define get_ndx(seq, do_store_orientation) ({\
    if ((seq.p & 1) || do_store_orientation > 1) {\
        fprintf(stderr, "orient already set at %s:%u\n",__FILE__,__LINE__);\
        fflush(stderr);\
        raise(SIGTRAP);\
    }\
    seq.t = seq.dna ^ seq.rc;\
    seq.t &= -seq.t;                  /* isolate deviant bit */\
    seq.t |= !seq.t;                  /* for palindromes: have to set one */\
    seq.t = !(seq.t & seq.dna);       /* was devbit set? */\
    seq.p |= do_store_orientation & seq.t;\
    seq.t = seq.dna ^ (-seq.t & (seq.dna ^ seq.rc)); /* dna or rc dependent on devbit */\
    seq.t ^= (-!!(seq.t & KEYNT_BUFSZ)) & SNDX_TRUNC_MASK; /*shorten index by one */\
})

static uint32_t
_build_ndx_kct(keyseq_t &seq, uint8_t const*const s, uint32_t with_orient = 1)
{
    uint32_t p = seq.p = _b2pos_of(seq.p);
    build_seq(s, seq, p);
    get_ndx(seq, with_orient);
    if (seq.t >= KEYNT_BUFSZ) {
        fprintf(stderr, "_build_ndx_kct(): seq.t(%u) >= KEYNT_BUFSZ\n", seq.t);
        fflush(stderr);
        raise(SIGTRAP);
    }
    return seq.t;
}

/* With the b6 conversion, only the specified characters are converted to
 * 2bits, left shifted by one, with all other bits zeroed. With the function
 * below one can test whether the ascii character prior to conversion was a
 * valid convertable character.
 */
#define isb6(b) (((b) | B6_MASK) == B6_MASK)

/* A similar, not entirely reversible conversion of ascii to b6. It ignores case
 * and works for both DNA and RNA, but all characters in the range [@-GP-W`-gp-w]
 * produce values interpreted as b6, including 'R', which occurs on hg19, chr3.
 */
# define __qb6(c) (c ^ ((c & (B6_ALT_CASE|B6_RNA)) | (B6_A2B ^ B6_B2A)))

#endif // RK_B6_H


