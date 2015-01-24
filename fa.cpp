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
    uint8_t b;
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

            assert(key < KEYNT_BUFSZ);
            if (kc->kpndx[key] == ~0u) {
                //fprintf(stderr, "%x\n", key);
                _buf_grow(kc->kp, 1ul);
                kc->kp[kc->kp_l].ct = 0;
                kc->kpndx[key] = kc->kp_l++;
            }
            kcs *kcp = kc->kp + kc->kpndx[key];

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
                        goto out; // FIXME: skip Y for now: it occurs twice. (will be overwritten)
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
int prepare_keys(kct_t* kc, unsigned readlength)
{
    unsigned k, itnr = 0;
    unsigned uq = kc->reg_l, tot = 0;
    for (k = 0; k != kc->kp_l; ++k) {
        kcs* kcp = &kc->kp[k];
        if (kcp->ct == 1) {
            _buf_grow(kc->reg, 1ul);
            kc->reg[kc->reg_l].pos = kcp->lastp;
            kc->reg[kc->reg_l].nr = itnr;
            kc->reg[kc->reg_l].type = 0;
            kc->reg_l++;
        }
        tot++;
    }
    uq = kc->reg_l - uq;
    fprintf (stdout, "%u / %u unique keys\n", uq, tot);
    return 0;
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
    kc.kp = _buf_init(kc.kp, 16);
    kc.kpndx = _buf_init_arr(kc.kpndx, KEYNT_BUFSZ_SHFT);
    for (i = 0; i != KEYNT_BUFSZ; ++i)
        kc.kpndx[i] = ~0u;

    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }

    res = fa_kc(&kc, g, gc);
    if (res < 0) { fprintf(stderr, "== index failed:%d", res);  goto out; }
    res = prepare_keys(&kc, seq->readlength);
    if (res < 0) { fprintf(stderr, "== counting failed:%d", res);  goto out; }

    ret = res != 0;

out:
    _buf_free(kc.kp);
    _buf_free(kc.kpndx);
    _buf_free(kc.reg);
    _buf_free(kc.hdr);
    _buf_free(kc.seq);
    kh_destroy(UQCT, kc.H);
    return ret;

}


