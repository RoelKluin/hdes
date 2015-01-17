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
inline uint64_t revseq(uint64_t x);
inline unsigned get_twisted_even(unsigned rev, unsigned sq);

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


