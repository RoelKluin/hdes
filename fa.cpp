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
    return fa->s != NULL ? 0 : -ENOMEM;
}

static void free_ndx(seqb2_t *fa)
{
    if (fa->s) free(fa->s);
    fa->s = NULL;
}

static inline int increment_keyct(uint8_t* p, uint64_t key, kct_t* kct)
{
    if (*p != 255) {
        (*p)++;
    } else {
        // the max keyct is sticky.
        // hash lookup and increment
        typeof(kct->H) H = kct->H;
        khiter_t k = kh_get(UQCT, H, key);
        if (k == kh_end(H)) {
            int err;
            k = kh_put(UQCT, H, key, &err);
            if_ever (err < 0) return err;
        }
        kh_val(H, k)++;
    }
    return 0;
}

static inline int
write_kcpos(uint8_t* p, uint64_t key, kct_t* kct)
{
    if (*p == 0) {
        fprintf(stderr, "0x%lx == %u, assertion fails at %u!\n", key, *p, kct->l);
        return -1;
    }
    kct->x[kct->l] = *p;

    if (++kct->l == BUF_STACK) {
	int c = gzwrite(kct->out, kct->x, kct->l);
	if (c < 0) return c - 1;
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

    while ((c = gc(g)) != '>') if (c == -1 || c == '@') return c; // skip to first
    do {
        while (c != '\n' && ((c = gc(g)) >= 0)) fputc(c, stderr); /* header etc */

        while (c != '>' && (c > 0)) { // foreach stretch of N's
            while((c = gc(g)) == 'N' || isspace(c)) {};

            // (re)initialize key
            register uint64_t t = 0, dna = 0ul, rev = 0ul;
            while (t != KEY_WIDTH && likely(c != '>' && c != -1)) {
                if (!isspace(c)) {
                    add_b6N0(c, c, dna, rev);
                    ++t;
                }
                c = gc(g);
            }

            // t is 2nd cNt bit excised, rev or dna
            t = get_b2cn_key(t, dna, rev);
            // a keycount of 255 == 0xff is a sticky max.
            t = kc->process(s + t, t, kc);

            while ((dna || c != 'N') && c != '>' && c != -1 && t == 0) {
                if (!isspace(c)) {
                    // seq or revcmp according to 2nd bit of central Nt
                    add_b6N0(t, c, dna, rev);
                    t = get_b2cn_key(t, dna, rev);
                    t = kc->process(s + t, t, kc);
                }
                c = gc(g);
            }
            if (t) c = (int)t;
        }
    } while (c > 0);

    return c & -(c != -1);
}

int
read_s(struct gzfh_t* fhin, struct seqb2_t* seq, kct_t* kc)
{
    int c;
    char *s = (char*)seq->s;
    uint64_t m = seq->m;
    fprintf(stderr, "==reading buffer s\n");
    while (m > INT_MAX) {
        c = gzread(fhin->io, s, INT_MAX);
        fprintf(stderr, "==read %d bytes\n", c);
        if (c == -1) return -1;
        m -= INT_MAX;
        s += INT_MAX;
    }
    c = gzread(fhin->io, s, m);
    fprintf(stderr, "==read %d bytes\n", c);
    if (c == -1) return -1;

    // kh_init already allocated struct.
    if (gzread(fhin->io, &kc->H->n_buckets, sizeof(khint_t)) == -1)
        return -1;

    if (gzread(fhin->io, &kc->H->size, sizeof(khint_t)) == -1)
        return -1;

    if (gzread(fhin->io, &kc->H->n_occupied, sizeof(khint_t)) == -1)
        return -1;

    if (gzread(fhin->io, &kc->H->upper_bound, sizeof(khint_t)) == -1)
        return -1;

    khint_t nb = kc->H->n_buckets;
    size_t t = __ac_fsize(nb) * sizeof(khint32_t);
    kc->H->flags = (khint32_t*)krealloc(kc->H->flags, t);
    if (kc->H->flags == NULL) return -ENOMEM;
    if (gzread(fhin->io, kc->H->flags, t) == -1)
        return -1;

    t = nb * sizeof(*kc->H->keys);
    kc->H->keys = (typeof(kc->H->keys))krealloc(kc->H->keys, t);
    if (kc->H->keys == NULL) return -ENOMEM;
    if (gzread(fhin->io, kc->H->keys, t) == -1)
        return -1;

    t = nb * sizeof(*kc->H->vals);
    kc->H->vals = (typeof(kc->H->vals))krealloc(kc->H->vals, t);
    if (kc->H->vals == NULL) return -ENOMEM;
    if (gzread(fhin->io, kc->H->vals, t) == -1)
        return -1;
    fprintf(stderr, "==done reading khash\n");

    return 0;
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
    kc.H = kh_init(UQCT);

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
            fhin2->name = &file[0]; // use as input 2 instead.
            fhin2->fd = fhout->fd;
            fhin2->fp = fhout->fp;

            res = set_io_fh(fhin2, blocksize, 0);
            fprintf(stderr, "== reopened %s:%d\n", fhin2->name, res);
            if (res < 0) break;

            res = read_s(fhin2, seq, &kc);
            fprintf(stderr, "== reread %s:%d\n", fhin2->name, res);
            if (res < 0) break;

            fprintf(stderr, "== closing %s\n", fhin2->name);
            if (fhin2->close && fhin2->close(fhin2->io) != Z_OK)
                fprintf(stderr, "main: gzclose fails for %s\n", fhin2->name);
            close(fhin2->fd);
            fhin2->fp = NULL; // prevent segfault on close in main()
            continue;
        } else if (i == 0) {
            fprintf(stderr, "== %s does not yet exist\n", fhout->name);
            memset(seq->s, '\0', seq->m);
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
        if (res == 0 && i == 1) {
            kc.x[kc.l] = '\0';
            res = gzwrite(fhout->io, line, ++kc.l);
            if (res > 0) res = 0;
        }
        if (res < 0) { fprintf(stderr, "== .index failed:%d", res);  break; }

        if (i == 0) {
            // the kcpos index function alreay wrote to disk, no read needed here.
            khint_t nb = kc.H->n_buckets;
            size_t nf = __ac_fsize(nb) * sizeof(khint32_t);
            size_t nk = nb * sizeof(*kc.H->keys);
            size_t nv = nb * sizeof(*kc.H->vals);
            size_t kisz = sizeof(khint_t);
            size_t new_m = seq->m + 4 * kisz + nf + nk + nv;

            char* s = (char*)realloc(seq->s, new_m);
            if (s == NULL) {
                ret = -ENOMEM;
                goto out;
            }
            seq->s = (uint8_t*)s;
            s += seq->m;
            seq->m = new_m;

            fprintf(stderr, "== appending khash (%u)\n", nb);
            memcpy(s, &kc.H->n_buckets, kisz); s += kisz;
            memcpy(s, &kc.H->size, kisz); s += kisz;
            memcpy(s, &kc.H->n_occupied, kisz); s += kisz;
            memcpy(s, &kc.H->upper_bound, kisz); s += kisz;
            memcpy(s, kc.H->flags, nf); s += nf;
            memcpy(s, kc.H->keys, nk); s += nk;
            memcpy(s, kc.H->vals, nv);
            assert(seq->m == ((uint8_t*)s - seq->s) + nv);

            fprintf(stderr, "== writing to disk\n"); fflush(NULL);
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
    kh_destroy(UQCT, kc.H);
    free_ndx(seq);
    return ret;

}

