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
//#include <glib.h>
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
static int fa_keyct(seqb2_t *seq, void* g, int (*gc) (void*))
{
    register unsigned i = 0;
    register int c;
    register uint8_t* p;
    uint8_t* s = seq->s;

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
            t = i + KEY_WIDTH;
            while (i != t && likely(c != '>' && c != -1)) {
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
            // a keycount of 255 == 0xff is a sticky max.
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

static int fa_kcbed(seqb2_t *seq, void* g, int (*gc) (void*))
{
    register int c;
    register unsigned i = 0;
    register uint8_t* p;
    uint8_t* s = seq->s;
    unsigned lastct = 0;
    int t;
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    gzFile out = fhout->io;
    uint8_t tid[256];
    char line[BUF_STACK];

    while ((c = gc(g)) != '>') if (c == -1 || c == '@') return c; // skip to first
    do {
        uint8_t* x = tid;
        while (!isspace(c) && ((*x++ = c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */
        *x = '\0';
        x = (uint8_t*)line;

        while (c != '\n' && ((c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */

        i = 0;
        while (c != '>' && c > 0) { // foreach stretch of N's
            while((t = (c = gc(g)) == 'N') || isspace(c)) i += t;

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

            register uint64_t b = (dna & KEYNT_STRAND) ? rev : dna;
            // excise out 2nd cNt bit
            b = ((b >> 1) & ~HALF_KEYNT_MASK) | (b & HALF_KEYNT_MASK);
            p = s + (b & KEYNT_TRUNC_MASK);
            lastct = *p;
            assert(lastct != 0);
            x = sprints(x, tid);
            *x++ = '\t';
            x += sprintu(x, i, '\t'); // begin pos

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
                        x += sprintu(x, i, '\t'); // end pos
                        x += sprintu(x, lastct, '\n'); // and counts
                        c = (char*)x - line;
                        if (BUF_STACK - c < 64) {
                            *x = '\0';
                            c = gzwrite(out, line, c);
                            if (c < 0) return c;
                            x = (uint8_t *)line;
                        }
                        lastct = *p;
                        x = sprints(x, tid); // next entry...
                        *x++ = '\t';
                        x += sprintu(x, i, '\t');
                    }
                }
                c = gc(g);
            }
            if (c == -1) c = 0;
            x += sprintu(x, i, '\t');
            x += sprintu(x, lastct, '\n');
            t = gzwrite(out, line, (char*)x - line);
            x = (uint8_t *)line;
            if (t < 0) c = t;
        }

        fprintf(stderr, "processed %u Nts\n", i + 1);
    } while (c > 0);

    return c;
}

int
read_s(seqb2_t *seq)
{
    struct gzfh_t* fhin2 = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    int c;
    char* s = (char*) seq->s;
    while ((c = gzread(fhin2->io, s, INT_MAX)) > 0) {
        if (c < 0) break;
        fprintf(stderr, "==%d bytes read\n", c);
        c /= sizeof(char), s += c;
        if (c != INT_MAX) return 0;
    }
    if (gzeof(fhin2->io) == 1) return -1;
    fprintf(stderr, "%s\n", gzerror(fhin2->io, &c));
    return -EIO;
}


int fa_print(seqb2_t *fa)
{
    return 1;
}

int
fn_convert(struct seqb2_t* seq, const char* search, const char* replace)
{
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    if (fhout->fp)
        fclose(fhout->fp);

    char* f = strstr(fhout->name, search);
    if (f == NULL) return 0;
    strncpy(f, replace, strlen(replace) + 1);
    return 1;
}

int fa_index(struct seqb2_t* seq)
{
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    struct gzfh_t* fhin = seq->fh + 2;
    uint64_t blocksize = (uint64_t)seq->blocksize << 20;
    int i, res, ret = -1;
    int is_gzfile = fhin->io != NULL;
    char file[256];
    void* g;
    int (*gc) (void*);
    int (*index_func) (seqb2_t*, void*, int (*) (void*));

    struct index_action ndxact[2] = {
        {".fa", ".keyct.gz", &fa_keyct},
        {".keyct.gz", ".kcbed.gz", &fa_kcbed}
    };

    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }

    if (fhout->name != NULL && ((res = strlen(fhout->name)) > 255)) {
        strncpy(file, fhout->name, res);
    } else {
        res = strlen(fhin->name); // create enough mem for keyct name
        if (strncmp(fhin->name, "stdin", res) == 0) {
            fputs("Cannot index from stdin (seek)\n", stderr);
            goto out;
        }
        strcpy(file, fhin->name);
    }
    fhout->name = &file[0];
    init_ndx(seq);

    for (i = 0; i != ARRAY_SIZE(ndxact) && res >= 0; ++i) {
        index_func = ndxact[i].index;

        if (!fn_convert(seq, ndxact[i].search, ndxact[i].replace))
            continue; // skip if this file was not requested on the commandline
        fprintf(stderr, "== %s\n", fhout->name);
        fhout->fp = fopen(fhout->name, "r"); // test whether file exists
        // XXX reading from hg19.keyct.gz file gives a different hg19.kcbed.gz than when preprocessed in mem XXX
        // hg19.keyct.gz files are reproducible
        if (fhout->fp) {
            res = set_io_fh(fhout, blocksize, 2); //2: force read
            if (res < 0) break;
            res = read_s(seq);
        } else {
            res = set_io_fh(fhout, blocksize, 0);
            if (res < 0) break;

            res = index_func(seq, g, gc);
            if (res < 0) break;

            if (is_gzfile)
                gzrewind(fhin->io);
            else
                rewind(fhin->fp);

            if (i == 0) {// XXX
                res = fhout->write(fhout, (char*)seq->s, ~0u, sizeof(char));
                if (fhout->close && fhout->close(fhout->io) != Z_OK)
                    fprintf(stderr, "main: gzclose fails for %s\n", fhout->name);
                close(fhout->fd);
            }
        }
    }
    ret = res != 0;

out:
    free_ndx(seq);
    return ret;

}

