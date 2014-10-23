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
#include <assert.h>
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

/*
 * Create key. Up to MULTIMAPPER just count. Thereafter insert/update hash.
 * hash contains an index for kct->mm. Each has 9 uints. 8 for each nucleotide
 * and strand orientation combination. This enables internal linkage, and one
 * count.
 */
static inline int
increment_keyct(seqb2_t *seq, kct_t* kct)
{
    // key is 2nd cNt bit excised, rev or dna

    int t = kct->dna & KEYNT_STRAND;
    unsigned key = t ? kct->dna : (kct->rev ^ 0xaaaaaaaa);
    key = (((key >> 1) & ~HALF_KEYNT_MASK) |
            (key & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;

    uint8_t* p = seq->s + key;
    khiter_t k;
    khash_t(UQCT) *H = kct->H;
    t = *p == _MULTIMAPPER;
    if (t == NO_MULTIMAPPER) {
        if (_MULTIMAPPER < 16) {
        // could instead use only 1, 2 or 4 bits only to account multimappers.
        // then 2nd middle Nt bit should be excised and determine nibble to incr.
        } else if (++(*p) == _MULTIMAPPER) {
            k = kh_get(UQCT, H, key);
            assert (k == kh_end(H)); // must be new.
            t = NEW_MULTIMAPPER;
        }
    }
    if (t != NO_MULTIMAPPER) {
        unsigned l = NO_LINK_YET;
        unsigned *n = &l;
        if (kct->last_mmpos + 1u == kct->pos) { // adjacent multimapper
           // handle internal linkage:
            l = ((kct->dna)& 3) | (-!(kct->dna & (KEYNT_STRAND << 2)) & 4);
            n = kct->mm + kct->last + l; // update linked instead
            l = *n;
        }

        if (l == NO_LINK_YET) {
            if (t == NEW_MULTIMAPPER) {
                int err;
                k = kh_put(UQCT, H, key, &err);
                if_ever (err < 0) return err;
                l = *n = kh_val(H, k) = kct->mm_l;
                kct->mm_l += MM_CT + 1;
                if (kct->mm_l >= kct->mm_m) {
                    kct->mm_m <<= 1;
                    unsigned* mm = (unsigned*)realloc(kct->mm,
                            kct->mm_m * sizeof(unsigned));
                    if (mm == NULL) return -ENOMEM;
                    kct->mm = mm;
                }
                for (int i = 0; i != MM_CT; ++i)
                    kct->mm[l + i] = NO_LINK_YET;
                kct->mm[l + MM_CT] = _MULTIMAPPER;
            } else {
                k = kh_get(UQCT, H, key);
                assert (k != kh_end(H));
                l = *n = kh_val(H, k);
            }
        } else assert(t == EXISTING_MULTIMAPPER);
        kct->last = l;
        kct->mm[l + MM_CT]++;
        kct->last_mmpos = kct->pos;
    }
    return 0;
}

static inline int
write_kcpos(seqb2_t *seq, kct_t* kct)
{
    int t = kct->dna & KEYNT_STRAND;
    unsigned key = t ? kct->dna : (kct->rev ^ 0xaaaaaaaa);
    key = (((key >> 1) & ~HALF_KEYNT_MASK) |
            (key & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;
    uint8_t* p = seq->s + key;
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    if (*p == 0) {
        fprintf(stderr, "ERROR: new key %x at %u in reiteration for kcpos!\n",
                key, kct->pos);
        return -1;
    }
    if (*p == _MULTIMAPPER) {
        unsigned l = NO_LINK_YET;
        khash_t(UQCT) *H = kct->H;
        unsigned *n = &l;
        t = kct->last_mmpos + 1u == kct->pos; // adjacent multimapper
        if (t) { // there may be some internal linkage left
            t = ((kct->dna)& 3) | (-!(kct->dna & (KEYNT_STRAND << 2)) & 4);
            n = kct->mm + kct->last + t; // update linked instead
            l = *n;
            // make t >= 0 for mmap ends: beyond seq->readlength - KEY_LENGTH
            t = kct->at[kct->at_l] - kct->pos - (seq->readlength - KEY_LENGTH);
        }
        if (t > 0) { // start of non-multimapper region, moves along.
            kct->at[kct->at_l] = kct->pos + 1;
        } else if (t == 0) { // start or end
            if (kct->at_l + 1 > kct->at_m) {
                unsigned* at = (unsigned*)realloc(kct->at,
                        kct->at_m * sizeof(unsigned));
                if (at == NULL) return -ENOMEM;
                kct->at = at;
            }
            if (kct->last_mmpos == kct->pos) {
                kct->at_l += 2; // initiate next start of non-multimapper region.
                kct->at[kct->at_l] = kct->pos + 1u;
            } else {
                kct->at[kct->at_l + !!kct->last_mmpos] = kct->pos;
            }
        }
        if (l == NO_LINK_YET) {
            khiter_t k = kh_get(UQCT, H, key);
            if (k == kh_end(H)) {
                fprintf(stderr, "k == kh_end(H): key:%x, at %u!\n", key, kct->pos);
                return -EINVAL;
            }
            l = *n = kh_val(H, k);
        }
        kct->last = l;
        kct->last_mmpos = kct->pos;
    }
    kct->x[kct->l] = *p;

    if (++kct->l == BUF_STACK) {
	int c = gzwrite(fhout->io, kct->x, kct->l);
	if (c < 0) return c - 1;
        kct->l = 0;
    }
    return 0;
}

/*
 * zero on success.
 */
static int fa_kc(seqb2_t *seq, kct_t* kc, void* g, int (*gc) (void*))
{
    register int c;
    register uint64_t b, t = KEY_WIDTH + 255;
    kc->dna = kc->rev = 0ul;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first
    c &= -(c != '>');

    while (c == 0) {
        c = gc(g);
        if ((t <= KEY_WIDTH) && !isspace(c)) {

            /* convert to twobit */
            b = c ^ ((c | B6_UC) & B6_LC);
            b ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);
            b ^= -((c & 0x8e) == 0x4) & B6_RNA;
            //fprintf(stderr, "%c:%lu\n", c, t);

            /* is twobit? => c zero*/
            c &= -((b | B6_MASK) != B6_MASK);
            if (c) { /* exceptions: header, N's, end or odd Nts */
                if (c == '>') {
                    kc->last_mmpos = 0u;
                    kc->pos = -1u; /* later incremented. */
                    t = KEY_WIDTH + 256; /* later decremented: start with header */
                } else {
                    if (c != 'N') {
                        if (c == -1) break;
                        fprintf(stderr, "Ignoring strange nucleotide:%c\n", c);
                    }
                    if (kc->dna == 0ul) /* rebuild key before next entry */
                        t = KEY_WIDTH;
                }
                b = c ^= c; // enter twobit 0 (A), and prevent error-out
            }
            /* append the twobit to dna and rev */
            kc->dna = (kc->dna << 2) | (b >> 1);
            kc->dna &= KEYNT_MASK; // prevent setting carriage bit.
            kc->rev = (b << (KEYNT_TOP - 1)) | (kc->rev >> 2);
            if (t == 0ul) { // is key complete?
                c = kc->process(seq, kc);
                if (c < 0)
                    fprintf(stderr, "process fails\n");
            }
            else --t;
            ++(kc->pos);
        } else {
            if (t > KEY_WIDTH) {
                t -= KEY_WIDTH;
                fputc(c, stderr);
                // '\0' for each space or last
                kc->tid[255 - t] = c & -!(isspace(c) | (t == 0));
                // iterate header up to newline, but fill only until 255th char 
                t = KEY_WIDTH + ((t - (t != 0)) & -(c != '\n'));
            }
            c ^= c;
        }
    }
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

    if (gzread(fhin->io, &kc->mm_l, sizeof(unsigned)) == -1)
        return -1;

    t = kc->mm_l * sizeof(unsigned);
    kc->mm = (unsigned*)krealloc(kc->mm, t);
    if (kc->mm == NULL) return -ENOMEM;
    if (gzread(fhin->io, kc->mm, t) == -1)
        return -1;
    fprintf(stderr, "==done reading mmap\n");

    return 0;
}


int fa_print(seqb2_t *fa)
{
    return 1;
}

int
fn_convert(struct gzfh_t* fhout, const char* search, const char* replace)
{
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
    kct_t kc = { .process = &increment_keyct};
    kc.H = kh_init(UQCT);
    kc.mm_m = 1 << 8;
    kc.mm_l = 0;
    kc.at_m = 1 << 8;
    kc.at_l = 0;

    unsigned* mm = (unsigned*)malloc(kc.mm_m * sizeof(unsigned));
    if (mm == NULL) return -ENOMEM;
    kc.mm = mm;

    unsigned* at = (unsigned*)malloc(kc.at_m * sizeof(unsigned));
    if (at == NULL) return -ENOMEM;
    kc.at = at;


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

        if (!fn_convert(fhout, ndxact[i], ndxact[i + 1])) //XXX still needed?
            continue; // skip if this file was not requested on the commandline

        fprintf(stderr, "== %s(%lu)\n", fhout->name, fhout - seq->fh);
        fhout->fp = fopen(fhout->name, "r"); // test whether file exists
        if (fhout->fp) {
            fprintf(stderr, "== %s already exists\n", fhout->name);
            if (i == 1) {
                // only keyct should be read in seq->s.
                fprintf(stderr, "== done.");
                break;
            }

            struct gzfh_t* fhin2 = seq->fh;
            assert(fhin2->name == NULL);
            fhin2->name = &file[0]; // use as 2nd input instead.
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
                fprintf(stderr, "gzclose fails for %s\n", fhin2->name);

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
            kc.l = 0;
        }
        res = fa_kc(seq, &kc, g, gc);
        if (res == 0 && i == 1) {
            kc.x[kc.l] = '\0';
            res = gzwrite(fhout->io, line, ++kc.l);
            if (res > 0) res = 0;
        }
        if (res < 0) { fprintf(stderr, "== index failed:%d", res);  break; }

        if (i == 0) {
            // the kcpos index function alreay wrote to disk, no read needed here.
            khint_t nb = kc.H->n_buckets;
            size_t nf = __ac_fsize(nb) * sizeof(khint32_t);
            size_t nk = nb * sizeof(*kc.H->keys);
            size_t nv = nb * sizeof(*kc.H->vals);
            size_t kisz = sizeof(khint_t);
            size_t usz = sizeof(unsigned);
            size_t n_mm = usz * kc.mm_l;
            size_t new_m = seq->m + 4 * kisz + nf + nk + nv + usz + n_mm;

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
            memcpy(s, kc.H->vals, nv); s += nv;
            fprintf(stderr, "== appending mmap (%lu)\n", n_mm);
            memcpy(s, &kc.mm_l, usz); s += usz;
            memcpy(s, kc.mm, n_mm);
            assert(seq->m == (uint8_t*)s - seq->s + n_mm);

            fprintf(stderr, "== writing to disk\n"); fflush(NULL);
            res = fhout->write(fhout, (char*)seq->s, seq->m, sizeof(char));
            if (res < 0) { fprintf(stderr, "== ->write failed:%d", res);  break; }

            fprintf(stderr, "== closing %s\n", fhout->name); // last before valg
            if (fhout->close && fhout->close(fhout->io) != Z_OK)
                fprintf(stderr, "gzclose fails for %s\n", fhout->name);
            close(fhout->fd);
            fhout->fp = NULL;

            if (is_gzfile) {// valgrind complaind here about uninit value.
                gzclearerr(fhin->io); // XXX: why needed???
                fprintf(stderr, "rewinding gz.\n");
                gzrewind(fhin->io);
            } else
                rewind(fhin->fp);
        }
    }
    ret = res != 0;

out:
    free(kc.mm);
    kh_destroy(UQCT, kc.H);
    free_ndx(seq);
    return ret;

}


