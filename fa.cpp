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
#include <stdlib.h> // realloc()
#include <ctype.h> // isspace()
#include <stdio.h> //fprintf()
#include <errno.h> // ENOMEM
#include <string.h>
#include <limits.h> //INT_MAX
#include "fa.h"

/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
static inline int
init_ndx(seqb2_t *fa)
{
    fa->m = sizeof(*fa->s) * KEYNT_BUFSZ;
    fa->s = (uint8_t*)malloc(fa->m);
    if (fa->s == NULL) return -ENOMEM;
    memset(fa->s, '\0', fa->m);

    return 0;
}

/*
 * zero on success.
 */
int fa_ndx(seqb2_t *fa)
{
    register int c;
    fputs("==Initializing memory", stderr);
    if ((c = init_ndx(fa)) < 0) {
        fprintf(stderr, "ERROR: init_ndx() returned %d\n", c);
        return c;
    }

    uint8_t* s = fa->s;
    void* g;
    int (*gc) (void*);
    if (fa->fh[2].io) {
        g = fa->fh[2].io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fa->fh[2].fp;
        gc = (int (*)(void*))&fgetc;
    }
    register unsigned i = 0;
    register uint8_t* p;

    while ((c = gc(g)) != '>') if (c == -1 || c == '@') return c; // skip to first
    do {
        while (c != '\n' && ((c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */

        i = 0;
        while (c != '>' && c > 0) { // foreach stretch of N's
            unsigned t = i;
            while((c = gc(g)) == 'N' || isspace(c)) i += (c == 'N');
            t = i - t;
            if (t) fprintf(stderr, "Skipped %d N's\n", t);

            // (re)initialize key
            register uint64_t dna = 0ul, rev = 0ul;
            while (i != KEY_WIDTH && likely(c != '>' && c != -1)) {
                if (!isspace(c)) {
                    c = b6(c);
                    c = -isb6(c) & (c >> 1); // zero [min] for N - but also for A.

                    dna = (dna << 2) | c;
                    rev = ((uint64_t)c << KEYNT_TOP) | (rev >> 2); // xor 0x2 after loop
                    ++i;
                }
                c = gc(g);
            }
            rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement

            register uint64_t b = (dna & KEYNT_STRAND) ? rev : dna;
            // excise out 2nd cNt bit
            b = ((b >> 1) & ~HALF_KEYNT_MASK) | (b & HALF_KEYNT_MASK);
            p = s + (b & KEYNT_TRUNC_MASK);
            if (*p != 255) (*p)++;

            while ((dna || c != 'N') && c != '>' && c != -1) {
                if (!isspace(c)) {
                    // seq or revcmp according to 2nd bit of central Nt
                    b = b6(c);
                    b = -isb6(b) & (b >> 1);

                    dna = (dna << 2) | b;
                    rev = ((b ^ 2) << KEYNT_TOP) | (rev >> 2);
                    b = (dna & KEYNT_STRAND) ? rev : dna;
                    b = ((b >> 1) & ~HALF_KEYNT_MASK) | (b & HALF_KEYNT_MASK);
                    p = s + (b & KEYNT_TRUNC_MASK);
                    if (*p != 255) (*p)++;
                    ++i;
                }
                c = gc(g);
            }
            if (c == -1) c = 0;
        }
        fprintf(stderr, "processed %u Nts\n", i + 1);
    } while (c > 0);
    if (c == 0) {
        struct gzfh_t* fh = &fa->fh[ARRAY_SIZE(fa->fh) - 1];
        size_t l = sizeof(*s) * ~0u;
        while (l && c >= 0) {
            c = l > INT_MAX ? INT_MAX : l;
            if ((c = fh->write(fh, (char*)s, c)) >= 0) l -= c;
            s += c * sizeof(char)/ sizeof(*s);
        }
        c = c < 0;
    }
    free(fa->s);
    fa->s = NULL;
    return c;
}

int fa_print(seqb2_t *fa)
{
    return 1;
}

