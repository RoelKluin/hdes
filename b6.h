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
#define KEY_CUTOFF 0                        // cut off, to get random keys, [XXX: keep for assembly]
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF) // entire key maximized on (after conversion)
#define KEYNT_TOP ((KEY_WIDTH - 1) * 2)

// maximum length of the read name
#define SEQ_MAX_NAME_ETC    (1u << 14)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ_SHFT ((KEY_LENGTH * 2) - 1)
#define KEYNT_BUFSZ (1ul << KEYNT_BUFSZ_SHFT)
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1ul) // also kctndx max.

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_MASK ((1u << (KEY_WIDTH << 1)) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)
#define KEYNT_TRUNC_UPPER (~HALF_KEYNT_MASK & KEYNT_TRUNC_MASK)


#define SNDX_TRUNC_MASK (KEYNT_TRUNC_MASK | KEYNT_BUFSZ)
#define HALF_SNDX_MASK (HALF_KEYNT_MASK | KEYNT_STRAND)
#define SNDX_TRUNC_UPPER (~HALF_SNDX_MASK & SNDX_TRUNC_MASK)



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
        dna = (dna << 2) | t;\
        rev = ((uint64_t)t << KEYNT_TOP) | (rev >> 2);\
})

#define add_b(t, dna, rev) ({\
        dna = (dna << 2) | t;\
        rev = ((uint64_t)t << KEYNT_TOP) | (rev >> 2);\
})

# define ischar_ignore_cs(c, c2) (((c ^ c2) & ~B6_ALT_CASE) == 0)

unsigned b6(unsigned c);
unsigned b6_spec(unsigned c, unsigned cs, unsigned no_u);
#define REVSEQ(dna, m) ({\
    m = 0x3333333333333333;\
    dna = ((dna & m) << 2) | ((dna & ~m) >> 2);\
    m = 0x0f0f0f0f0f0f0f0f;\
    dna = ((dna & m) << 4) | ((dna & ~m) >> 4);\
    m = 0x00ff00ff00ff00ff;\
    asm ("bswap %0" : "=r" (dna) : "0" (dna));\
    dna >> (64 - (KEY_WIDTH << 1));\
})

inline uint64_t revcmp(uint64_t dna) /* Reverse Complement, is ok. */
{
    uint64_t m;
    dna ^= 0xaaaaaaaaaaaaaaaaUL;
    return REVSEQ(dna, m);
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


