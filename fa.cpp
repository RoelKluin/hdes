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

static int
insert_mmapper(kct_t* kct, unsigned* l, uint64_t key)
{
    int t;
    khiter_t k = kh_put(UQCT, kct->H, key, &t);
    if_ever (t < 0) return t;
    kh_val(kct->H, k) = *l = kct->mm_l;
    _buf_grow(kct->mm, 1);
    kct->mm_l ++;
    return 0;
}

/*
 * process fasta, store 2bit string and account keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*))
{
    uint32_t dna = 0u, rev = 0u, key, pos = 0u;
    int c;
    uint16_t t = KEY_WIDTH;
    uint8_t b2, b;
    uint8_t *s = kc->seq;
    *s = '\0';

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID
    if_ever(gzungetc(c, (gzFile)g) == -1) return -EINVAL;

    while (1) {
        c = gc(g);
        if ((t <= KEY_WIDTH) && !isspace(c)) {

            /* convert to twobit */
            b = c ^ ((c | B6_UC) & B6_LC);
            b ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);
            b ^= -((c & 0x8e) == 0x4) & B6_RNA;

            /* is twobit? => zero c and skip character check branch */
            c &= -((b | B6_MASK) != B6_MASK);
            if (c) { /* the exceptions are: N's or odd Nts, headers and EOF */
                if (c == '>') {
                    t = REF_CHANGE;
                    _buf_grow(kc->reg, 1ul);
                    kc->reg[kc->reg_l].pos = pos;
                    kc->reg[kc->reg_l].nr = kc->hdr_l;
                    kc->reg[kc->reg_l].type = 2;
                    continue;
                }
                if (c == -1)
                    return 0;
                _buf_grow(kc->reg, 1ul);
                kc->reg[kc->reg_l].pos = pos;
                kc->reg[kc->reg_l].type = 1;
                key = 0;
                t = N_STRETCH;

                if (c != 'N')
                    fprintf(stderr, "Ignoring strange nucleotide:%c\n", c);

                continue;
            }
            b >>= 1;
            dna = (dna << 2) & KEYNT_MASK; // do not set carriage bit.
            dna |= b;
            rev = (b << KEYNT_TOP) | (rev >> 2);

            *s |= b;
            if ((++pos & 3) == 0) {
                _buf_grow2(kc->seq, 1ul, s);
                *++s = '\0';
                kc->seq_l++;
            }
            *s <<= 2;
            if (t) {
                --t;
                continue;
            }
            // key is 2nd cNt bit excised, rev or dna
            key = (dna & KEYNT_STRAND) ? dna : (rev ^ 0xaaaaaaaa);
            key = (((key >> 1) & KEYNT_TRUNC_UPPER) | (key & HALF_KEYNT_MASK));

            kcs* kcp = &kc->kp[key];
            kcp->lastp = pos; // replace with new l
            kcp->ct++;
        } else {
            if (t == N_STRETCH) {
                if (c != 'N') {
                    if (c == '\n') continue;
                    kc->reg[kc->reg_l].nr = key;
                    kc->reg_l++;
                    if (c == -1) return 0;
                    if_ever(gzungetc(c, (gzFile)g) == -1) return -EINVAL;
                    t = KEY_WIDTH;
                }
                ++key;
            } else if (t == REF_CHANGE) {
                fputc(c, stderr);
                // '\0' for each space or last
                _buf_grow(kc->hdr, 1ul);
                kc->hdr[kc->hdr_l++] = !isspace(c) ? c : '\0';
                // iterate header up to newline, but fill only until 255th char
                if (c == '\n') {
                    if (kc->reg_l && kc->hdr[kc->reg[kc->reg_l].nr] == 'Y')
                        goto out; // FIXME: skip Y for now: it occurs twice.
                    kc->reg_l++;
                    t = KEY_WIDTH;
                }
            }
        }
    }
out:
    kc->seq_l = pos;
    return 0;
}


// process keys and occurance in each range
// FIXME: currently only prints each
void prepare_keys(kct_t* kc)
{
    unsigned k;
    char s[KEY_LENGTH+1];
    for (k=0; k != KEY_LENGTH; ++k) s[k] = 'A';
    s[k] = '\0';
    k = 0;
    fprintf (stdout, "kseq\t\tlastpos\t\tcount\tkeynr\n");
    while (1) {
        kcs* kcp = &kc->kp[k];

        // next key: 1. calculate parity
        unsigned c;
        if (__builtin_parity(k) & 1) {
            c = (k & -k) << 1;
            if (c >= KEYNT_BUFSZ) break; // wrap
        } else {
            c = 1u;
        }
        fprintf (stdout, "%s\t%u\t%u\t0x%x\n", s, kcp->lastp, kcp->ct, k);
        k ^= c;
        c = __builtin_ctz(c);
        c >>= 1;
        s[c] = b6(((k >> (c << 1)) & 3) << 1);
    }
}

int
fa_print(seqb2_t *fa)
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
    //char line[BUF_STACK];

    const char* ndxact[4] = {".fa", "_kcpos.wig.gz"};

    if (fhout->name != NULL && ((res = strlen(fhout->name)) > 255)) {
        strncpy(file, fhout->name, res);
    } else {
        res = strlen(fhin->name);
        if (strncmp(fhin->name, "stdin", res) == 0) {
            fputs("Cannot index from stdin (seek)\n", stderr);
            return -1;
        }
        strcpy(file, fhin->name);
    }
    fhout->name = &file[0];

    kct_t kc;
    kc.H = kh_init(UQCT);
    kc.seq = _buf_init(kc.seq, 16);
    kc.hdr = _buf_init(kc.hdr, 8);
    kc.reg = _buf_init(kc.reg, 8);
    kc.kp = _buf_init(kc.kp, KEYNT_BUFSZ_SHFT);


    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }
    memset(kc.kp, 0, sizeof(kc.kp[0]) << kc.kp_m);

    res = fa_kc(&kc, g, gc);
    if (res < 0) { fprintf(stderr, "== index failed:%d", res);  goto out; }
    prepare_keys(&kc);

    ret = res != 0;

out:
    _buf_free(kc.kp);
    _buf_free(kc.reg);
    _buf_free(kc.hdr);
    _buf_free(kc.seq);
    kh_destroy(UQCT, kc.H);
    return ret;

}


