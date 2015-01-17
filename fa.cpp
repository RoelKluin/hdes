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
static int
write_kcwig_fixedStep0(seqb2_t *seq, kct_t* kct)
{
    if (kct->l + 512 >= BUF_STACK) {
        struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
	int c = gzwrite(fhout->io, kct->x, kct->l);
	if (c < 0) return c;
        kct->l = 0;
    }
    uint8_t* w = (uint8_t*)&kct->x[kct->l];
    w = sprints(w, (uint8_t*)"fixedStep\tchrom=chr");
    w = sprints(w, (uint8_t*)kct->tid);
    // wig is one-based, each key spans KEY_WIDTH Nt's, first ia from 1-KEYWIDTH
    w = sprints(w, (uint8_t*)"\tstart=1");

    // FIXME: this is fasta specific, but needed
    // value between 4th and 5th ':' in tid is start
    uint8_t* s = (uint8_t*)kct->tid;
    for(int i=0; *++s != ':' || ++i != 4;) {}
    while (*++s != ':') *w++ = *s;

    // span must be 1 or wigtobigwig fails. Alternatively write bed file
    w = sprints(w, (uint8_t*)"\tstep=1\tspan=1\n");
    //w += sprintu(w, KEY_WIDTH, '\n');
    kct->l += w - (uint8_t*)&kct->x[kct->l];
    return 0;
}

static int
write_kcpos(seqb2_t *seq, kct_t* kct)
{
    uint64_t key = (kct->dna & KEYNT_STRAND) ? kct->dna :
        (kct->rev ^ 0xaaaaaaaaaaaaaaaa);
    key = (((key >> 1) & KEYNT_TRUNC_UPPER) | (key & HALF_KEYNT_MASK));

    unsigned l = seq->s[key];
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    if (l == 0 && kct->Nmask == 0u) {
        fprintf(stderr, "ERROR: new key %lx at %u in reiteration for kcpos!\n",
                key, kct->pos);
        return -1;
    }
    if (l == _MULTIMAPPER || (kct->Nmask && kct->last_mmpos + 1u == kct->pos)) {
        l = UNDEFINED_LINK;
        int t = kct->last_mmpos + 1u == kct->pos; // adjacent multimapper
        if (t) { // there may be some internal linkage left
            l = kct->last + (((kct->dna)& 3) | (-!(kct->dna & (KEYNT_STRAND << 2)) & 4));
            l = kct->mm[l];
            // make t >= 0 for mmap ends
            t = kct->at[kct->at_l] - kct->pos;
        }
        if (t > 0) { // start of non-multimapper region, moves along.
            kct->at[kct->at_l] = kct->pos + 1;
        } else if (t == 0) { // start or end
            _buf_grow(kct->at, 1);
            if (kct->last_mmpos == kct->pos) {
                kct->at_l += 2; // initiate next start of non-multimapper region.
                kct->at[kct->at_l] = kct->pos + 1u;
            } else {
                // kct->last_mmpos == 0: is start of chromosome.
                kct->at[kct->at_l + !!kct->last_mmpos] = kct->pos;
            }
        }
        if (l == UNDEFINED_LINK) {
            khiter_t k = kh_get(UQCT, kct->H, key);
            if (k == kh_end(kct->H)) {
                if (kct->Nmask) {
                    assert (kct->last_mmpos + 1u == kct->pos);
                    t = insert_mmapper(kct, &l, key);
                    if_ever (t < 0) return t;
                    kct->mm[l + MM_CT] = seq->s[key];
                } else {
                    fprintf(stderr, "k == kh_end(H): key:%lx, at %u!\n", key, kct->pos);
                    return -EINVAL;
                }
            }
            l = kh_val(kct->H, k);
        }
        kct->last = l;
        l = kct->mm[l + MM_CT];
        kct->last_mmpos = kct->pos;
    }
    if (kct->l + 16 >= BUF_STACK) {
	int c = gzwrite(fhout->io, kct->x, kct->l);
	if (c < 0) return c - 1;
        kct->l = 0;
    }
    if (kct->header == &write_kcwig_fixedStep0) // FIXME
        kct->l += sprintu((uint8_t*)&kct->x[kct->l], l, '\n');
    else
        kct->x[kct->l] = seq->s[key];
    return 0;
}*/

/*
 * zero on success.
 */
/*static int
fa_kc2(seqb2_t *seq, kct_t* kc, void* g, int (*gc) (void*))
{
    uint64_t dna = 0ul, rev = 0ul;
    uint16_t t = KEY_WIDTH;
    uint8_t b2, b;

    for (unsigned i = 0; i != kc->pos; ++i) {
        if ((i & 3) == 0)
            b2 = kc->seq[kc->pos >> 2];
        b = b2 & 3;
        b2 >>= 2;

        // append the twobit to dna and rev
        kc->dna = (kc->dna << 2) | b;
        kc->dna &= KEYNT_MASK; // prevent setting carriage bit.
        kc->rev = (b << KEYNT_TOP) | (kc->rev >> 2);
        if (t == 0ul) { // is key complete?
            c = kc->process(seq, kc);
            if (c < 0)
               fprintf(stderr, "process fails\n");
        } else {
            --t;
        }
    }
    return c & -(c != -1);
}*/

/*
 * process fasta, store 2bit string and account keys.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*))
{
    uint32_t dna = 0u, rev = 0u, key, dl, pos;
    int c;
    uint16_t t = 0xffff, Nmask = 0u;
    uint8_t b2, b;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID
    if (c == '>') c = 0;

    while (1) {
        c = gc(g);
        if ((t <= KEY_WIDTH) && !isspace(c)) {

            /* convert to twobit */
            b = c ^ ((c | B6_UC) & B6_LC);
            b ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);
            b ^= -((c & 0x8e) == 0x4) & B6_RNA;

            Nmask = (Nmask << 1) & N_MASK;
            /* is twobit? => zero c and skip character check branch */
            c &= -((b | B6_MASK) != B6_MASK);
            if (c) { /* the exceptions are: N's or odd Nts, headers and EOF */
                if (c == '>') {
                    Nmask = 0u;
                    t = 0xffff;
                    _buf_grow(kc->ho, 1ul);
                    kc->ho[kc->ho_l++] = kc->hdr_l;
                    continue;
                }
                if (c == -1) return 0;
                Nmask |= 1u;
                if (Nmask == N_MASK && (pos & 3) == 0) /* only N's? rebuild key */
                    t = KEY_WIDTH + 2; // 2: + 1 but is decremented later..

                if (c != 'N')
                    fprintf(stderr, "Ignoring strange nucleotide:%c\n", c);

                b = '\0';
            }
            b >>= 1;
            dna = (dna << 2) | b;
            dna &= KEYNT_MASK; // prevent setting carriage bit.
            rev = (b << KEYNT_TOP) | (rev >> 2);

            b = b << ((~pos & 3) << 1);
            if ((pos & 3) == 0) {
                if (pos) {
                    _buf_grow(kc->seq, 1ul);
                    kc->seq[kc->seq_l++] = b2;
                }
                b2 = b;
            }
            b2 |= b;
            ++pos;
            if (t == 0 && Nmask == 0u) {
                // key is 2nd cNt bit excised, rev or dna
                key = (dna & KEYNT_STRAND) ? dna : (rev ^ 0xaaaaaaaa);
                key = (((key >> 1) & KEYNT_TRUNC_UPPER) | (key & HALF_KEYNT_MASK));

                kcs* kcp = &kc->kp[key];
                dl = kc->seq_l - kcp->lastp; // distance
                // no bits required to store distance
                dl = (sizeof(dl) << 3) - __builtin_clz(dl);
                if (dl > kcp->dbits) kcp->dbits = dl;

                kcp->lastp = kc->seq_l; // replace with new l

                // account for several distances the counts
                if (dl <= 8) {
                    if(kcp->ct1 == 0xfffffff) {
                        fprintf(stderr, "0-8: count reached max for key %x\n", key);
                        return -EINVAL;
                    }
                    kcp->ct1++;
                    // TODO: may want to skip 2bit storage for such repeats
                } else if (dl <= 16) {
                    if(kcp->ct2 == 0xfffff) {
                        fprintf(stderr, "8-16: count reached max for key %x\n", key);
                        return -EINVAL;
                    }
                    kcp->ct2++;
                } else if (dl <= 24) {
                    if(kcp->ct3 == 0xffff) {
                        fprintf(stderr, "16-24: count reached max for key %x\n", key);
                        return -EINVAL;
                    }
                    kcp->ct3++;
                } else {
                    if(kcp->ct4 == 0xffffff) {
                        fprintf(stderr, "24-32: count reached max for key %x\n", key);
                        return -EINVAL;
                    }
                    kcp->ct4++;
                }
            }
            t -= (t != 0);
        } else {
            if (Nmask == N_MASK) {
                if (c != 'N') {
                    if (c == '\n') continue;
                    if (c == -1) return 0;
                    _buf_grow(kc->ho, 1ul);
                    kc->ho[kc->ho_l++] = pos | 0x80000000u;
                    t = KEY_WIDTH;
                    if_ever(gzungetc(c, (gzFile)g) == -1) return -EINVAL;
                } else {
                    ++pos;
                }
            } else if (t > KEY_WIDTH) {
                fputc(c, stderr);
                // '\0' for each space or last
                _buf_grow(kc->hdr, 1ul);
                kc->hdr[kc->hdr_l++] = !isspace(c) ? c : '\0';
                // iterate header up to newline, but fill only until 255th char
                if (c == '\n') {
                    if (kc->ho_l && kc->hdr[kc->ho[kc->ho_l - 1]] == 'Y')
                        return 0; // FIXME: skip Y for now: it occurs twice.
                    t = KEY_WIDTH;
                    pos = 0u;
                    dna = 0ul;
                    rev = 0ul;
                } else {
                    --t;
                    assert(t != KEY_WIDTH);
                }
            }
        }
    }
    return 0;
}

/*
procedure add (n: integer; A,B:word; PA,PB:bit;
               var S:word; var PS:bit; var CE, CF:bit);
var i: integer; E, F, T: bit;
begin
   E := PA; F := PB;
   for i:= 0 to n-1 do begin {in parallel, using previous inputs}
       S[i] := (E and F) ⊕ A[i] ⊕ B[i];
       E := (E and (not F)) ⊕ A[i];
       F := ((not E) and F) ⊕ B[i];
   end;
   CE := E; CF := F;
end;
 */

// process keys and occurance in each range
// FIXME: currently only prints each
void prepare_keys(kct_t* kc)
{
    unsigned k;
    char s[KEY_LENGTH+1];
    for (k=0; k != KEY_LENGTH; ++k) s[k] = 'A';
    s[k] = '\0';
    k = 0;
    fprintf (stdout, "kseq\t\tlastpos\t\tmaxdbit\t<=8\t<=16\t<=24\t<=32"
            "\ttot\tkeynr\n");
    while (1) {
        kcs* kcp = &kc->kp[k];
        unsigned tot = (unsigned)kcp->ct1 + (unsigned)kcp->ct2 +
            (unsigned)kcp->ct3 + (unsigned)kcp->ct4;

        // next key: 1. calculate parity
        unsigned c;
        if (__builtin_parity(k) & 1) {
            c = (k & -k) << 1;
            if (c >= KEYNT_BUFSZ) break; // wrap
        } else {
            c = 1u;
        }
        fprintf (stdout, "%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t0x%x\n", s,
                (unsigned)kcp->lastp, (unsigned)kcp->dbits, (unsigned)kcp->ct1,
                (unsigned)kcp->ct2, (unsigned)kcp->ct3, (unsigned)kcp->ct4, tot, k);
        k ^= c;
        c = __builtin_ctz(c);
        c >>= 1;
        s[c] = b6(((k >> (c << 1)) & 3) << 1);
    }
}

/*
 * Create key. Up to MULTIMAPPER just count. Thereafter insert/update hash.
 * hash contains an index for kct->mm. Each has 9 uints. 8 for each nucleotide
 * and strand orientation combination. This enables internal linkage, and one
 * count.
 */
/*static int
fa_kc_bad(seqb2_t *seq, kct_t* kc, void* g, int (*gc) (void*))
{
    uint64_t b, t = KEY_WIDTH + 255;
    uint64_t dna = 0ul, rev = 0ul, key;
    uint32_t Nmask;
    int c;
    unsigned pos, l;
    uint8_t b2, shift, bits;
    khiter_t k;
    unsigned __int128 locs;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first
    c &= -(c != '>');

    while (c == 0) {
        c = gc(g);
        if ((t <= KEY_WIDTH) && !isspace(c)) {

            // convert to twobit
            b = c ^ ((c | B6_UC) & B6_LC);
            b ^= (((c ^ 'U') & ~B6_ALT_CASE) != 0);
            b ^= -((c & 0x8e) == 0x4) & B6_RNA;

            // is twobit? => zero c and skip character check branch
            c &= -((b | B6_MASK) != B6_MASK);
            if (c) { // the exceptions are: header, N's, end or odd Nts
                if (c == '>') {
                    last_mmpos = Nmask = last = 0u;
                    pos = -1u; // later incremented.
                    t = KEY_WIDTH + 256; // later decremented: start with header
                    _buf_grow(kc->ho, 2ul);
                    kc->ho[kc->ho_l++] = pos + 1;
                    kc->ho[kc->ho_l++] = kc->hdr_l;
                } else {
                    Nmask |= 1u;
                    if (c != 'N') {
                        if (c == -1) return 0;
                        fprintf(stderr, "Ignoring strange nucleotide:%c\n", c);
                    }
                }
                b = c ^= c; // enter twobit 0 (A), and no need to error-out
            }

            // append the twobit to dna and rev
            b >>= 1;
            dna = (dna << 2) | b;
            dna &= KEYNT_MASK; // prevent setting carriage bit.
            rev = (b << KEYNT_TOP) | (rev >> 2);

            if (t) goto store_2bit; // skip if key cannot yet be complete

            if (Nmask)
                goto store_2bit; // if key uncertain, don't add it

            // key is 2nd cNt bit excised, rev or dna
            key = (dna & KEYNT_STRAND) ? dna : (rev ^ 0xaaaaaaaaaaaaaaaa);
            key = (((key >> 1) & ~HALF_KEYNT_MASK) |
                    (key & HALF_KEYNT_MASK)) & KEYNT_TRUNC_MASK;
            locs = seq->s[key];

            if ((locs & ULL(1)) == ULL(0)) { // set when multimapper
                // try to add diminuitive position in locs.
                p = 0u;
                // TODO: partial insertion?
                for (shift = 1, bits = 32; shift + bits < 127;) {
                    t = (locs & ((ULL(1) << (bits + shift)) - ULL(1))) >> shift;
                    if (t == 0) {
                        seq->s[key] |= ULL((pos - p) << shift);
                        goto store_2bit; // fits, don't store in hash
                    } else if ((t > 16) || (p > )){
                        shift += bits;
                        p += t;
                        // set bits required for remaining positions
                        bits = (GSZ_MAX - p) - 1u;
                        bits |= bits >> 1;
                        bits |= bits >> 2;
                        bits |= bits >> 4;
                        bits |= bits >> 8;
                        bits |= bits >> 16;
                        ++bits;
                    } else { // repeat, next 8 bits is repeatcount
                        
                    }

                }
                (1ul << bits) - 1ul
            }
            if (l != _MULTIMAPPER && ++(seq->s[key]) != _MULTIMAPPER)
                goto store_2bit;

            k = kh_get(UQCT, kc->H, key);

            if (k == kh_end(kc->H)) { // new multimapper
                k = kh_put(UQCT, kct->H, key, &c);
                if_ever(c < 0) return c;
                c ^= c;

                kh_val(kct->H, k) = l = kc->mm_l;
                _buf_grow(kct->mm, 1);
                kc->mm[kct->mm_l++] = _MULTIMAPPER; // we encountered this many already
            } else {
                l = kh_val(kc->H, k);
            }
            kc->mm[l]++;
store_2bit:
            t -= !!t; // key is one Nt more complete now (if not allready).
            b <<= (~pos & 3) << 1;
            if ((pos & 3) == 0) {
                if (pos) {
                    _buf_grow(kc->seq, 1ul); // XXX
                    kc->seq[kc->seq_l++] = b2;
                }
                b2 = b;
            }
            b2 |= b;
            ++pos;
            Nmask = (Nmask << 1) & N_MASK;
        } else {
            if (t > KEY_WIDTH) {
                t -= KEY_WIDTH;
                fputc(c, stderr);
                // '\0' for each space or last
                _buf_grow(kc->hdr, 1ul);
                hdr[hdr_l++] = c & -!(isspace(c) | (t == 0));
                // iterate header up to newline, but fill only until 255th char
                t = KEY_WIDTH + ((t - (t != 0)) & -(c != '\n'));
                if (t == KEY_WIDTH && kc->hdr[kc->ho[kc->ho_l - 1]] == 'Y')
                    return 0; // FIXME: skip Y for now: it occurs twice.
            }
            c ^= c;
        }
    }
    return c & -(c != -1);
}*/

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
    kc.mm = _buf_init(kc.mm, 8);
    kc.at = _buf_init(kc.at, 8);
    kc.seq = _buf_init(kc.seq, 16);
    kc.hdr = _buf_init(kc.hdr, 8);
    kc.ho = _buf_init(kc.ho, 2);
    kc.kp = _buf_init(kc.kp, KEYNT_BUFSZ_SHFT);


    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }
    assert(sizeof(kc.kp[0]) == 16);
    memset(kc.kp, 0, 16ul << kc.kp_m);

    res = fa_kc(&kc, g, gc);
    if (res < 0) { fprintf(stderr, "== index failed:%d", res);  goto out; }
    prepare_keys(&kc);

/*FIXME:
    if (!fn_convert(fhout, ndxact[0], ndxact[1]))
        goto out; // skip if this file was not requested on the commandline

    fprintf(stderr, "== %s(%lu)\n", fhout->name, fhout - seq->fh);
    fhout->fp = fopen(fhout->name, "r"); // test whether file exists
    if (fhout->fp) {
        fprintf(stderr, "== %s already exists\n== done.\n", fhout->name);
        goto out;
    }

    res = set_io_fh(fhout, blocksize, 0);
    if (res < 0) { fprintf(stderr, "== set_io_fh failed:%d", res); goto out; }

// FIXME:
    kc.process = &write_kcpos;
    kc.header = &write_kcwig_fixedStep0; // XXX remove to write binary instead
    kc.x = line;
    kc.l = 0;

    res = fa_kc2(seq, &kc, g, gc);
// XXX
    if (res == 0) {
        kc.x[kc.l] = '\0';
        res = gzwrite(fhout->io, line, ++kc.l);
        if (res > 0) res = 0;
    }
    if (res < 0) { fprintf(stderr, "== index failed:%d", res);  goto out; }
*/
    ret = res != 0;

out:
    _buf_free(kc.kp);
    _buf_free(kc.ho);
    _buf_free(kc.hdr);
    _buf_free(kc.seq);
    _buf_free(kc.at);
    _buf_free(kc.mm);
    kh_destroy(UQCT, kc.H);
    return ret;

}


