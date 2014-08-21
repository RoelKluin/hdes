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
#include <stdio.h> // fprintf(), fseek();
#include <errno.h> // ENOMEM
#include <string.h> // memset()
#include <limits.h> //INT_MAX
#include <sys/types.h>
#include <unistd.h>
#include "fa.h"

/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
static int
init_ndx(seqb2_t *fa)
{
    fa->m = sizeof(*fa->s) * KEYNT_BUFSZ;
    fa->s = (uint8_t*)malloc(fa->m);
    if (fa->s == NULL) return -ENOMEM;
    memset(fa->s, '\0', fa->m);

    return 0;
}

static void free_ndx(seqb2_t *fa)
{
    if (fa->s) free(fa->s);
    fa->s = NULL;
}

/*
 * zero on success.
 */
static int fa_ndx(uint8_t* s, void* g, int (*gc) (void*))
{
    register int c;
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
            t = 0;
            while (t != KEY_WIDTH && likely(c != '>' && c != -1)) {
                if (!isspace(c)) {
                    c = b6(c);
                    c = -isb6(c) & (c >> 1); // zero [min] for N - but also for A.

                    dna = (dna << 2) | c;
                    rev = ((uint64_t)c << KEYNT_TOP) | (rev >> 2); // xor 0x2 after loop
                    ++t;
                }
                c = gc(g);
            }
            i += t;
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

    return c;
}

static int fa_bed(uint8_t* s, void* g, int (*gc) (void*), int fno)
{
    register int c;
    register unsigned i = 0;
    register uint8_t* p;
    char line[256];
    unsigned lastct = 0;
    int t;

    while ((c = gc(g)) != '>') if (c == -1 || c == '@') return c; // skip to first
    do {
        char* en, *st = &line[0];
        en = st;
        while (!isspace(c) && ((*st++ = c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */
        *st++ = '\t';

        while (c != '\n' && ((c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */

        i = 0;
        while (c != '>' && c > 0) { // foreach stretch of N's
            while((c = gc(g)) == 'N' || isspace(c)) i += (c == 'N');

            // (re)initialize key
            register uint64_t dna = 0ul, rev = 0ul;
            t = 0;
            while (t != KEY_WIDTH && likely(c != '>' && c != -1)) {
                if (!isspace(c)) {
                    c = b6(c);
                    c = -isb6(c) & (c >> 1); // zero [min] for N - but also for A.

                    dna = (dna << 2) | c;
                    rev = ((uint64_t)c << KEYNT_TOP) | (rev >> 2); // xor 0x2 after loop
                    ++t;
                }
                c = gc(g);
            }
            rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement
            i += t;
            en = st + sprintu((uint8_t *)st, i, '\t');

            register uint64_t b = (dna & KEYNT_STRAND) ? rev : dna;
            // excise out 2nd cNt bit
            b = ((b >> 1) & ~HALF_KEYNT_MASK) | (b & HALF_KEYNT_MASK);
            p = s + (b & KEYNT_TRUNC_MASK);
            lastct = *p;

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
                    ++i;
                    if (*p != lastct) {
                        en += sprintu((uint8_t *)en, i, '\t');
                        en += sprintu((uint8_t *)en, lastct, '\n');
                        c = write(fno, line, en - line);
                        if (c < 0) return -EINVAL;
                        en = st + sprintu((uint8_t *)st, i, '\t');
                    }
                }
                c = gc(g);
            }
            if (c == -1) c = 0;
        }
        en += sprintu((uint8_t *)en, i, '\t');
        en += sprintu((uint8_t *)en, lastct, '\n');
        t = write(fno, line, en - line);
        if (t < 0) c = t;

        fprintf(stderr, "processed %u Nts\n", i + 1);
    } while (c > 0);

    return c;
}

int
rd_keyct(seqb2_t *seq)
{
    struct gzfh_t* fhin2 = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    return 1;
}


int fa_print(seqb2_t *fa)
{
    return 1;
}

int
reopen(struct seqb2_t* seq, const char* search, const char* replace, int (*rd) (struct seqb2_t*))
{
    uint64_t blocksize = (uint64_t)seq->blocksize << 20;
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    if (fhout->fp)
        fclose(fhout->fp);

    char* f = strstr(fhout->name, search);
    if (f == NULL) return -EINVAL;
    strcpy(f, replace);
    fhout->fp = fopen(fhout->name, "r"); // test whether file exists
    if (fhout->fp) {
        fprintf(stderr, "Read from %s\n", fhout->name);
        if (rd)
            return rd(seq); // XXX XXX: read in preexisting data
        else
            return -EINVAL;
    } else {
        fprintf(stderr, "Write to %s\n", fhout->name);
        return set_io_fh(fhout, blocksize, 0);
    }
}

int fa_index(struct seqb2_t* seq)
{
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    struct gzfh_t* fhin = seq->fh + 2;
    int res, ret = -1, use_refname = fhout->name == NULL;
    int is_gzfile = fhin->io != NULL;

    if (use_refname) {
        res = strlen(fhin->name); // create enough mem for keyct name
        if (strncmp(fhin->name, "stdin", res) == 0) {
            fputs("Cannot index from stdin (seek)\n", stderr);
            goto out;
        }
        char* f = (char*)malloc(res + 20);
        if (f == NULL) {
            ret = -ENOMEM;
            goto out;
        }
        fhout->name = f;
        strncpy(f, fhin->name, res);
    }
    res = reopen(seq, ".fa", ".keyct.gz", &rd_keyct);
    if (res < 0) goto out;

    fputs("==Initializing memory", stderr);
    if (init_ndx(seq) < 0) {
        fputs("ERROR: init_ndx() failed\n", stderr);
        goto out;
    }

    fputs("Creating counts per key\n", stderr);
    if (is_gzfile) {
        res = fa_ndx(seq->s, fhin->io, (int (*)(void*))&gzgetc);
        if (res < 0) goto out;

        res = fhout->write(fhout, (char*)seq->s, ~0u, sizeof(char));
        if (res != 0) goto out;

        res = reopen(seq, ".keyct.gz", ".bed.gz", NULL);
        if (res != 0) goto out;

        gzrewind(fhin->io);
        res = fa_bed(seq->s, fhin->io, (int (*)(void*))&gzgetc, fileno(fhout->fp));
    } else {
        res = fa_ndx(seq->s, fhin->fp, (int (*)(void*))&fgetc);
        if (res < 0) goto out;

        res = fhout->write(fhout, (char*)seq->s, ~0u, sizeof(char));
        if (res != 0) goto out;

        res = reopen(seq, ".keyct.gz", ".bed.gz", NULL);
        if (res != 0) goto out;

        rewind(fhin->fp);
        res = fa_bed(seq->s, fhin->fp, (int (*)(void*))&fgetc, fileno(fhout->fp));
    }
    if (res == 0) ret = 0;;

out:
    free_ndx(seq);
    if (use_refname && fhout->name != NULL)
        free(fhout->name);
    return ret;

}

