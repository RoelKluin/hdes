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

static inline int increment_keyct(uint8_t* p, kct_t* nop)
{
    if (*p != 255) (*p)++; // the max keyct is sticky.
    return 0;
}

static inline int
write_kcpos(uint8_t* p, kct_t* kct)
{
    assert(*p != 0);
    kct->x[kct->l] = *p;

    if (++kct->l == BUF_STACK) {
	int c = gzwrite(kct->out, kct->x, kct->l);
	if (c < 0) return c;
        kct->l = 0;
    }
    return 0;
}

/*
 * zero or -1 on success.
 */
static int fa_kc(seqb2_t *seq, kct_t* kc, void* g, int (*gc) (void*))
{
    register int c;
    uint8_t* s = seq->s;
    int t;

    while ((c = gc(g)) != '>') if (c == -1 || c == '@') return c; // skip to first
    do {
        while (c != '\n' && ((c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */

        while (c != '>' && (c > 0)) { // foreach stretch of N's
            while((c = gc(g)) == 'N' || isspace(c)) {};

            // (re)initialize key
            register uint64_t b, dna = 0ul, rev = 0ul;
            t = 0;
            while (t != KEY_WIDTH && likely(c != '>' && c != -1)) {
                if (!isspace(c)) {
                    b = b6N0(c);
                    dna = (dna << 2) | b;
                    rev = ((uint64_t)b << KEYNT_TOP) | (rev >> 2); // xor 0x2 after loop
                    ++t;
                }
                c = gc(g);
            }
            rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement

            // excise out 2nd cNt bit
            b = get_b2cn_key(b, dna, rev);
            // a keycount of 255 == 0xff is a sticky max.
            t = kc->process(s + b, kc);
            if (t < 0) return t;

            while ((dna || c != 'N') && c != '>' && c != -1) {
                if (!isspace(c)) {
                    // seq or revcmp according to 2nd bit of central Nt
                    b = b6N0(c);
                    dna = (dna << 2) | b;
                    rev = ((b ^ 2) << KEYNT_TOP) | (rev >> 2);
                    b = get_b2cn_key(b, dna, rev);
                    t = kc->process(s + b, kc);
                    if (t < 0) return t;
                }
                c = gc(g);
            }
        }
    } while (c > 0);

    return c;
}

int
read_s(struct gzfh_t* fhin, char*s)
{
    int c;
    unsigned m = INT_MAX;
    while ((c = gzread(fhin->io, s, m)) == (int)m) {
        fprintf(stderr, "==%d bytes read\n", c);
        c /= sizeof(char);
        s += c;
    }
    fprintf(stderr, "==%d bytes read\n", c);
    if (c >= 0) return 0;
    if (gzeof(fhin->io) == 1) return -1;
    m = c;
    fprintf(stderr, "== err... %s:%d:%s\n", fhin->name, (int)m, gzerror(fhin->io, &c));
    fprintf(stderr, "== err:%d\n", c);
    return -EIO;
}


int fa_print(seqb2_t *fa)
{
    return 1;
}

int
fn_convert(struct gzfh_t* fhout, const char* search, const char* replace)
{
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
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    uint64_t blocksize = (uint64_t)seq->blocksize << 20;
    int i, res, ret = -1;
    void* g;
    int (*gc) (void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];
    char line[BUF_STACK];

    const char* ndxact[4] = {".fa", ".keyct.gz", ".kcpos.gz"};
    kct_t kc = { .process = &increment_keyct };


    if (fhout->name != NULL && ((res = strlen(fhout->name)) > 255)) {
        strncpy(file, fhout->name, res);
    } else {
        res = strlen(fhin->name);
        if (strncmp(fhin->name, "stdin", res) == 0) {
            fputs("Cannot index from stdin (seek)\n", stderr);
            goto out;
        }
        strcpy(file, fhin->name);
    }
    fhout->name = &file[0];
    init_ndx(seq);

    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }

    for (i = 0; i != 2 && res >= 0; ++i) {

        if (!fn_convert(fhout, ndxact[i], ndxact[i + 1]))
            continue; // skip if this file was not requested on the commandline
        fprintf(stderr, "== %s(%lu)\n", fhout->name, fhout - seq->fh);
        fhout->fp = fopen(fhout->name, "r"); // test whether file exists
        if (fhout->fp) {
            fprintf(stderr, "== %s already exists\n", fhout->name);
            if (i == 2 - 1) {
                // only keyct should be read in seq->s.
                fprintf(stderr, "== done.");
                break;
            }

            struct gzfh_t* fhin2 = seq->fh;
            assert(fhin2->name == NULL);
            fhin2->name = &file[0]; // us as input 2 instead.
            fhin2->fd = fhout->fd;
            fhin2->fp = fhout->fp;

            res = set_io_fh(fhin2, blocksize, 0);
            fprintf(stderr, "== reopened %s:%d\n", fhin2->name, res);
            if (res < 0) break;

            res = read_s(fhin2, (char*) seq->s);
            fprintf(stderr, "== reread %s:%d\n", fhin2->name, res);
            if (res < 0) break;

            fprintf(stderr, "== closing %s\n", fhin2->name);
            if (fhin2->close && fhin2->close(fhin2->io) != Z_OK)
                fprintf(stderr, "main: gzclose fails for %s\n", fhin2->name);
            close(fhin2->fd);
            fhin2->fp = NULL; // prevent segfault on close in main()
            continue;
        }

        res = set_io_fh(fhout, blocksize, 0);
        if (res < 0) { fprintf(stderr, "== set_io_fh failed:%d", res); break; }

        if (i == 1) {
            kc.process = &write_kcpos;
            kc.x = line;
            kc.out = fhout->io;
            kc.l = 0;
        }
        res = fa_kc(seq, &kc, g, gc);
        if (res == -1) {
            res = 0;
            if (i == 1) {
                kc.x[kc.l] = '\0';
                res = gzwrite(fhout->io, line, ++kc.l);
                if (res > 0) res = 0;
            }
        }
        if (res < 0) { fprintf(stderr, "== .index failed:%d", res);  break; }

        if (i == 0) {
            // the kcpos index function alreay wrote to disk, no read needed here.
            res = fhout->write(fhout, (char*)seq->s, seq->m, sizeof(char));
            if (res < 0) { fprintf(stderr, "== ->write failed:%d", res);  break; }
            fprintf(stderr, "== closing %s\n", fhout->name);

            if (fhout->close && fhout->close(fhout->io) != Z_OK)
                fprintf(stderr, "main: gzclose fails for %s\n", fhout->name);
            close(fhout->fd);
            fhout->fp = NULL;

            if (is_gzfile)
                gzrewind(fhin->io);
            else
                rewind(fhin->fp);
        }
    }
    ret = res != 0;

out:
    free_ndx(seq);
    return ret;

}

