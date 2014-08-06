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
//#include <stdio.h> //fprintf()
#include <errno.h> // ENOMEM
#include "fq.h"

/*
 * store seqphred and return corresponding 2bit Nt
 */
static inline unsigned
seqphred(uint8_t *b, int q)
{
    unsigned c = b6(*b);
    *b = q;

    // 1) seqphred includes redundant key Nts. We could save a few bytes per sequence.
    // 2) Maybe reads with N's should be processed at the end - in assembly.
    if (isb6(c)) {
        *b = (*b | (c << 5)) + 3;
        return c >>= 1;
    }
    return c ^ c; // zero [min] for N - but also for A.
}

static inline void
decode_key(uint8_t* seq, unsigned dna)
{
    uint8_t* revcmp = seq + (KEY_LENGTH << 1) + 1;
    *revcmp = '\0';
    //put back cNT bit in seqrc (zero)
    dna = ((dna & ~HALF_KEYNT_MASK) << 1) | (dna & HALF_KEYNT_MASK);

    unsigned i, c;
    for (i = 0; i != KEY_LENGTH; ++i, ++seq, dna >>= 2) {
        c = (dna & 0x3) << 1;
        *seq = b6(c ^ 4);
        *--revcmp = b6(c);
    }
    *seq = '|';
    assert(revcmp == seq + 1);
}

/* the following should be enough for 32 bit int */

static inline const uint8_t *
sprints(uint8_t *out, const uint8_t *s)
{
    while (*s) { *++out = *s; ++s; }
    return out;
}


static inline unsigned
sprintu(uint8_t *out, register uint64_t u)
{
    register uint8_t *s = out + 21;
    *s = '\0';

    while (u) {
        register unsigned t = u % 10;
        u /= 10;
        *--s = t + '0';
    }
    for (;*s; ++u, ++s) *++out = *s;
    return u;
}

static inline unsigned
sprint0x(uint8_t *out, register uint64_t u)
{
    *++out = '0'; *++out = 'x';
    register uint8_t *s = out + 16;
    *s = '\0';

    while (u) {
        register int t = u & 0xf;
        *--s = t + (t >= 10 ? 'a' - 10 : '0');
        u >>= 4;
    }
    for (u = 2;*s; ++u, ++s) *++out = *s;
    return u;
}

/**
 * decode the sequence and phreds and print the fastq entries.
 */
int
fq_print(seqb2_t *fq)
{
    const unsigned readlength = fq->readlength;
//    fprintf(stderr, "%u\n", readlength); fflush(NULL);
    const unsigned phred_offset = fq->phred_offset;
    uint8_t seqrc[(KEY_LENGTH+1)<<1];
    uint8_t *const s = fq->s + SEQ_OFFSET_BYTES + BUF_OFFSET_BYTES;

    register uint32_t j = KEYNT_BUFSZ - 1u;
    uint32_t* look = &fq->lookup[j];
    const unsigned bufm = 1u << 20;
    char buf[bufm];
    uint8_t *w = (uint8_t *)buf;
    int ret = 0;
    register unsigned c = 0;
    struct gzfh_t* fh = &fq->fh[ARRAY_SIZE(fq->fh) - 1];

    do {
        register uint32_t i = *look--;
        if (i != 1u) { // has keys
            decode_key(seqrc, j);
            do {
                c = w - (uint8_t *)buf;
                if ((c + (readlength << 1) + FQ_MAX_NAME_ETC) >= bufm) {
                    if (fh->write(fh, buf, c) < 0) {
                        ret = -1;
                        break;
                    }
                    w = (uint8_t *)buf;
                }

                register uint8_t *b = s + i - 1;
                i = *b; i <<= 8;
                i |= *--b; i <<= 8;
                i |= *--b; i <<= 8;
                i |= *--b; // next seq offset for same key, or 1u.
                // similarly key offset could be retrieved, but KEY_WIDTH
                // needs to be added to get to the real offset.
                b += BUF_OFFSET_BYTES;
                *w = '@'; *++w = *b;                     // write header
                while ((c = *++b)) *++w = c;
                *++w = ' ';  w += sprint0x(w, j);        // key in comment
                *++w = '\t'; w += sprints(w, seqrc) - w; // and regexp
                *++w = '\n';

                register uint8_t *d = w;
                w += readlength + 1; // jump to second header and write the '+' line
                *++w = '+'; *++w = '\n'; //XXX valgrind

                while((c = *++b) != 0xff) { // sequence and qual are printed here
                    if (c >= 3) {
                        c -= 3;
                        *++d = b6((c >> 5) & 0x6);
                        *++w = (c & 0x3f) + phred_offset;
                    } else {
                        *++d = 'N';
                        *++w = c + phred_offset;
                    }
                }
                *++d = *++w = '\n'; ++w; // for seq and qual
            } while (i != 1u);
        }
    } while (j-- != 0u);
    c = w - (uint8_t *)buf;

    if (c && fh->write(fh, buf, c) < 0) ret = -1;

    return ret;
}

/*
 * This function parses fastq lines
 *
 * Return value:
 *      > 0     wrong file format
 *      0       read entire file
 *     -1       truncated file
 *     -ENOMEM  allocation error
 */
int
fq_b2(seqb2_t *fq)
{
    uint64_t l = fq->l, m = fq->m;
    uint8_t *s = fq->s;
    uint32_t *const lookup = fq->lookup;
    unsigned nr = fq->nr, key_ct = fq->key_ct;
    const unsigned phred_offset = fq->phred_offset;
    unsigned fq_ent_max = (fq->readlength << 1) + FQ_MAX_NAME_ETC;
    register uint8_t *b = s + l;
    register int c;

    void* g;
    int (*gc) (void*);
    if (fq->fh[0].io) {
        g = fq->fh[0].io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fq->fh[0].fp;
        gc = (int (*)(void*))&fgetc;
    }

    while ((c = gc(g)) != '@') /* skip to first header */
        if (c == -1 || c == '>') goto out;
    do {
        register unsigned i = 0;
        assert((uint64_t)l + fq_ent_max < (1ul << 32));
        if (l + fq_ent_max >= m) { // grow buffer if insufficient space for another read
            m <<= 1;
            fprintf(stderr, "==realloc at %u reads to 0x%lx\n", nr, m);
            if ((s = (uint8_t *)realloc(fq->s, m)) == NULL) {
                c = -ENOMEM;
                goto out;
            }
            fq->s = s;
            b = s + l;
        }
        uint8_t* o = b;
        b += BUF_OFFSET_BYTES + SEQ_OFFSET_BYTES;
        while (!isspace(c = gc(g)) && (c >= 0)) { /* header */
            *b++ = c;
            if_ever(++i == FQ_MAX_NAME_ETC) { c = 2; goto out; }
        }
        while (c != '\n' && (c >= 0)) c = gc(g); /* comment */
        *b++ = '\0';
//        fprintf(stderr, "%s...\n", b);fflush(NULL);

        for (i = 0; (c = gc(g)) != '\n'; ++b, ++i) {
            if (c < 0) goto out;
            *b = c;
        }
        if_ever(i > fq->readlength) fq->readlength = i;
        *b = '\0';
        c = '+';
        if (gc(g) != c) goto out; // EOF

        while ((c = gc(g)) != '\n') // skip 2nd header
            if (c < 0) goto out;

//        fprintf(stderr, "'%s'\n", b); fflush(NULL);

//        fprintf(stderr, "%u, %u\n", i, KEY_WIDTH); fflush(NULL);
        if (i < KEY_WIDTH) { // sequence too short. skipped and qual as well
            b = o;
            while ((c = gc(g)) != '\n' && c != -1) --i;
            if (i) { c = -1; goto out; }
            c = gc(g);
            if (c >= 0) continue;
            break;
        }
        b -= i;
        i = KEY_WIDTH;
        if_ever (nr++ == fq->readlimit) break;

        // initialize key
        register uint64_t dna = 0ul, rev = 0ul;
        c = gc(g);

        do {
//            fprintf(stderr, "q:'%c'\n", c); fflush(NULL);
            if_ever ((c -= phred_offset) > 50) goto out;
            c = seqphred(b, c);
//            fprintf(stderr, "s:'%c'\n", b6(c<<1)); fflush(NULL);

            dna = (dna << 2) | c;
            rev = ((uint64_t)c << KEYNT_TOP) | (rev >> 2); // Note: xor 0x2 after loop
            c = gc(g); // get next quality
            ++b;
            if_ever (c == -1 || c == '\n') goto out;
        } while (--i != 0u);
        rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement

        // seq or revcmp according to 2nd bit of central Nt
        uint64_t lmax = (dna & KEYNT_STRAND) ? rev : dna;
        uint64_t maxv = lmax ^ (lmax >> 2); // maximize on gray Nt sequence
        unsigned maxi = i = 0;

        for (; *b != '\0'; ++i, ++b) {
//            fprintf(stderr, "2q:'%c'\n", c); fflush(NULL);
            if_ever ((c -= phred_offset) > 50) goto out;
            c = seqphred(b, c);
//            fprintf(stderr, "2s:'%c'\n", b6(c<<1)); fflush(NULL);

            // update window for both strands
            dna = ((dna << 2) & KEYNT_MASK & 0xffffffffffffffff) | c;
            rev = ((uint64_t)(c ^ 2) << KEYNT_TOP) | (rev >> 2);

            // select a strand dependent on the central Nt
            uint64_t key = (dna & KEYNT_STRAND) ? rev : dna;

            // with gray Nt code, we may select a more sequence specific section.
            uint64_t Ntgc = key ^ (key >> 2);

            if ((Ntgc > maxv) || unlikely(Ntgc == maxv && (dna & KEYNT_STRAND))) {
                // a new maximum, keep it.
                lmax = key;
                maxv = Ntgc;
                maxi = i;
            }

            if_ever ((c = gc(g)) == -1) {
                if (*b == '\0') c = 0; // ok: last fastq line without newline
                goto out; // get next phred score
            }
        }
//        fprintf(stderr, "2q:'%c'\n", c); fflush(NULL);
        if_ever ((c -= phred_offset) > 50) goto out;
        i = seqphred(b, c);
        *b = 0xff;
        assert(i == b6('A'));
        assert ((lmax & KEYNT_STRAND) == 0);
        lmax = ((lmax >> 1) & ~HALF_KEYNT_MASK) | (lmax & HALF_KEYNT_MASK); // excise out 2nd cNt bit
        lmax &= KEYNT_TRUNC_MASK; /* truncate to make keys more random */

        // FIXME: reorder SEQ_OFFSET and BUF_OFFSET
        for (i = 0; i != SEQ_OFFSET_BYTES; ++i, maxi >>= 8) // insert key offset
            *o++ = maxi & 0xff;

        uint32_t* r = lookup + lmax;
        i = *r;

        if (i == 1u) ++key_ct; // first element
        *r = l; // exchange by this offset

        for (c = 0; c != BUF_OFFSET_BYTES; ++c, i >>= 8)
            *o++ = i & 0xff; // move former buffer offset, or last element mark (1u) here

        l = ++b - s;
        c = gc(g);
    } while (c >= 0);
    c = 0;
out: fq->l = l, fq->m = m, fq->nr = nr, fq->key_ct = key_ct;
    return c;
}


