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
#include "fa.h"
#include "b6.h"

void
free_kc(kct_t *kc)
{
    buf_free(kc->bnd);
    buf_free(kc->id);
    buf_free(kc->s);
    buf_free(kc->kct);
    buf_free(kc->h);
    buf_free(kc->hkoffs);
    buf_free(kc->ext_iter);
    buf_free(kc->contxt_idx);
}

static void
print_kct(kct_t *kc, Mantra* at, Bnd &b, uint32_t* tk)
{
    unsigned i = 0;
    EPR("");
    Hdr* h = kc->h;
    uint8_t *s = kc->s;
    Mantra* bnd = kc->bnd_l ? kc->bnd : at;
    uint32_t* hkoffs = kc->hkoffs;

    for (uint32_t*k = kc->kct; k - kc->kct < kc->kct_l;) {
        if (bnd != b.obnd + b.obnd_l) {
            while (h - kc->h != bnd->ho)
                s += h++->len;
            hkoffs = kc->hkoffs + bnd->ho;
        } else {
            if (h - kc->h != kc->h_l - 1) {
                s += h++->len;
            } else {
                ++i; // increase extension nr.
                h = kc->h;
                s = kc->s;
            }
            if (++hkoffs == kc->hkoffs + kc->hkoffs_l)
                hkoffs = &kc->kct_l;
        }

        while (k - kc->kct < *hkoffs) {
            char c = k!=tk?(k!=b.tgtk?' ':'t'):'k';
            EPR0("%c>%s:%u\t%u\t%c", c, kc->id + h->ido, i, *hkoffs, *k & DUP_BIT?'*':' ');
            if (*k != 0)
                print_posseq(s, *k);
            else
                EPR("(removed)");

            ++k;
        }
        if (bnd != b.obnd + b.obnd_l && ++bnd == kc->bnd + kc->bnd_l)
            bnd = at;
    }
}

// if kct grows, pointers *thisk, k, b.tgtk become invalid
static void
buf_grow_ks(kct_t *kc, Bnd &b, unsigned add, uint32_t **k)
{
    if ((kc->kct_l + add) >= (1ul << kc->kct_m)) {
        EPR("buf_grow_ks grew");
        unsigned ok1, ok2;

        ok1 = *k - kc->kct;
        ok2 = b.tgtk - kc->kct;
        uint32_t *t = (uint32_t *)realloc(kc->kct, sizeof(uint32_t) << ++kc->kct_m);
        if_ever (t == NULL)
            raise(SIGTRAP);
        kc->kct = t;
        *k = t + ok1;
        b.tgtk = t + ok2;
    }
}

/*
 * b.obnd contains ranges per contig for sequence not yet be mappable, but may become so.
 * Initially there are one or more dependent on the presence and number of N-stretches;
 * initialized in key_init.cpp, But as adjacent uniques cover regions, the remaining `mantra'
 * shrinks (here).
 *
 * mantras start at the chr start, the position after a N-stretch or after a region scoped
 * by uniques, the first non-unique after. The mantra ends at chr end, stretch or first non
 * unique before a region scoped by uniques.
 *
 * while the no. regions can grow, the total no. sequence decreases, hence: shrink.
 */


static void
k_compression(kct_t *kc, Bnd &b, uint32_t *hkoffs, uint32_t *k)
{

    NB(hkoffs <= kc->hkoffs + kc->h_l);
    while (hkoffs != kc->hkoffs + kc->h_l)
        *hkoffs++ = b.tgtk - kc->kct;

    Hdr* h = kc->h;
    b.s = kc->s;

    while (k - kc->kct != kc->kct_l) {
        NB(hkoffs < kc->hkoffs + kc->hkoffs_l);
        while (k - kc->kct == *hkoffs) {
            if (h - kc->h != kc->h_l - 1) {
                b.s += h++->len;
            } else {
                h = kc->h;
                b.s = kc->s;
            }
            *hkoffs++ = b.tgtk - kc->kct;
        }
        if (*k) {
            if (k != b.tgtk) {
                keyseq_t seq = {.p = b2pos_of(*k) };
                kc->contxt_idx[build_ndx_kct(kc, seq, b.s, 0)] = b.tgtk - kc->kct;
                *b.tgtk = *k;
                *k ^= *k;//
            }
            ++b.tgtk;
        }
        ++k;
    }
    kc->kct_l = b.tgtk - kc->kct;
    while (hkoffs != kc->hkoffs + kc->hkoffs_l)
        *hkoffs++ = kc->kct_l;
}

// keys between two uniqs.
static inline void
excise_one(kct_t *kc, Bnd &b, uint32_t *k)
{
    keyseq_t seq = {0};
    NB(k < kc->kct + kc->kct_l);
    NB(k != kc->kct + kc->kct_l);
    //*k &= ~DUP_BIT; // or keys that were moved after extension still have their dup bit set.
    kc->kct[kc->kct_l] = seq.p = *k;
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s, 0);
    *k ^= *k;//
    *contxt_idx = kc->kct_l++;
    ++b.moved;//S;
}

static void
excise(kct_t *kc, Bnd &b, uint32_t *thisk)
{
    for (uint32_t *k = b.tgtk + b.moved; k < thisk; ++k)
        excise_one(kc, b, k);
}

//keys not in scope of unique. hot
static void
move_uniq(kct_t *kc, Bnd &b, C unsigned start, C unsigned pend)
{
    unsigned dna = 0, rc = 0, t, key_complete = start - 2, p = start - NT_WIDTH;
    do { // build key
        t = (b.s[p>>3] >> (p&6)) & 3;
        rc = ((rc << 2) & KEYNT_MASK) | (t ^ 2);
        dna = t << KEYNT_TOP | dna >> 2;
        p += 2;
    } while (p != key_complete);

    do {
        t = (b.s[p>>3] >> (p&6)) & 3;
        rc = ((rc << 2) & KEYNT_MASK) | (t ^ 2);
        dna = t << KEYNT_TOP | dna >> 2;

        // build complement independent index
        unsigned dev = dna ^ rc;
        t = dev & -dev;           /* isolate deviant bit */
        unsigned ori = !((t | !t) & dna); /* for palindromes: have to set one. was devbit set? */
        t = dna ^ (-ori & dev); /* dna or rc dependent on devbit */
        dev = t & KEYNT_BUFSZ;
        t ^= (-!!(t & KEYNT_BUFSZ)) & SNDX_TRUNC_MASK; /*shorten index by one */

        p += 2;
        uint32_t *k = kc->kct + kc->contxt_idx[t];
        if (k > b.tgtk) {
            //O; 1st occurance
            ++kc->ct;    // unique or decremented later.
            *k = p | ori;      // set new pos and strand, unset dupbit

            *b.tgtk = *k;
            *k ^= *k;//
            kc->contxt_idx[t] = b.tgtk - kc->kct;
            if (b.tgtk[b.moved])
                --b.moved;

            ++b.tgtk;
        } else if (k < b.tgtk) {
            //O; second+ occurance (may still be in scope)
            // after &&: do not set dupbit if previous key was within read scope (a cornercase)
            if ((~*k & DUP_BIT) && (k - kc->kct < b.fk || b2pos_of(*k) + b.ext < p)) {
                //P; no dup after all

                *k |= DUP_BIT;
                --kc->ct;
            }
        } else {

            //O; also 1st occurance
            ++kc->ct;    // unique or decremented later.
            *k = p | ori;      // set new pos and strand, unset dupbit

            ++b.tgtk;
        }
    } while (p != pend);
}

/*
 * Called for extensions 1+, primary uniques were already determined upon reading the fasta.
 * secondary uniques arise when all but one key are already covered by primary unique keys.
 *
 * 1) iterate over as of yet ambiguous keys, determine which ones have become unique.
 * 2) determine uniques and shrink remaining ambiguous regions (mantra) accordingly.
 *
 * keys are kept ordered. Ambiguous keys are kept ordered upon first occurance. unique keys are
 * isolated and ordered on extension, contig and position.
 *
 */
static void
ext_uq_iter(kct_t *kc, Bnd &b)
{
    uint32_t *k = kc->kct;
    uint32_t *hkoffs = kc->hkoffs;
    Mantra *bnd = kc->bnd;
    Hdr *h = kc->h;
    unsigned skctl = kc->kct_l;

    kc->bnd = b.obnd; // buffer will be reused
    b.obnd = bnd; // swap b.obnd and kc->bnd.
    b.obnd_l = kc->bnd_l;
    b.obnd_m = kc->bnd_m;

    unsigned t = b.obnd_l;
    kroundup32(t);
    kc->bnd_m = __builtin_ctz(t) + 1;
    kc->bnd_l = 0;

    buf_realloc(kc->bnd, kc->bnd_m);
    //bnd = b.obnd; (already)

    b.tgtk = k;
    b.s = kc->s;
    b.fk = k - kc->kct;
    b.moved = 0;//B; initial state

    do {
        NB(h - kc->h <= bnd->ho);
        while (h - kc->h != bnd->ho) {
            //~ also update header
            *hkoffs++ = b.tgtk - kc->kct;
            t = hkoffs - kc->hkoffs;
            buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
            hkoffs = kc->hkoffs + t;
            b.s += h++->len;
            b.fk = k - kc->kct;
        }
        NB(hkoffs == kc->hkoffs + bnd->ho);

        //k - kc->kct <= *hkoffs: may be untrue after key excision.
        //B; next region

        while ((k - kc->kct) < *hkoffs && b2pos_of(*k) < bnd->e) {

            if (IS_UQ(k)) { //~ uniq
                unsigned end = b2pos_of(*k);
                if (bnd->s + b.ext >= end) {
                    //P; in scope of start; excision
                    bnd->s = end + 2;
                    ++k;
                    buf_grow_ks(kc, b, (k - b.tgtk) - b.moved, &k);
                    excise(kc, b, k);
                    --k;
                } else {
                    //P; out of scope.
                    move_uniq(kc, b, bnd->s, end - 2);
                    if (end + 2 == bnd->e) { //k; prevent insertion & removal, + b.ext?
                        bnd->e = end;
                        ++k;
                        buf_grow_ks(kc, b, (k - b.tgtk) - b.moved, &k);
                        excise(kc, b, k);
                        break;
                    }
                    Mantra copy = *bnd;
                    copy.e = end;
                    buf_grow_add(kc->bnd, 1ul, 0, copy);
                    bnd->s = end + 2;
                    buf_grow_ks(kc, b, 1, &k);
                    excise_one(kc, b, k);
                }
            }
            ++k;
        }

        if (bnd->s + b.ext >= bnd->e) {
            //P; in scope of start; excision (no bnd copy from b.obnd to kc->bnd)
            buf_grow_ks(kc, b, (k - b.tgtk) - b.moved, &k);
            excise(kc, b, k);
        } else if (b2pos_of(*k) != bnd->e) { //XXX: hg19: assertion '*k != 0' failed
            //P; alt
            move_uniq(kc, b, bnd->s, bnd->e - 2);
            buf_grow_add(kc->bnd, 1ul, 0, *bnd);//K;
        }
    } while (++bnd != b.obnd + b.obnd_l);

    t = hkoffs - kc->hkoffs;
    while (h - kc->h != kc->h_l) {
        buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
        ++h;
    }
    k_compression(kc, b, kc->hkoffs + t, k);//K;

    NB(skctl == kc->kct_l, "skctl(%u) != kc->kct_l(%u)", skctl, kc->kct_l);//B; final state
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res;
    Bnd b = {0};
    unsigned end = (kc->readlength - KEY_WIDTH + 1) << 1;
    kc->ext_iter = buf_init(kc->ext_iter, 1);
    for (b.ext = 0; b.ext != end; b.ext += 2) {
        kc->iter = 0;
        if (kc->hkoffs[kc->h_l-1] != 0) { // or all keys were already finished.
            do { // until no no more new uniques
                kc->ct = 0;
                ext_uq_iter(kc, b);
                EPR("observed %u potential in iteration %u, extension %u\n",
                    kc->ct, ++kc->iter, b.ext >> 1);
            } while (kc->ct > 0);
        }
        EPR("----[ end of extension %u ]-------", b.ext >> 1);
        buf_grow_add(kc->ext_iter, 1ul, 0, kc->iter);
    }
    _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
    _ACTION(save_kc(fhout + 3, kc), "writing unique keycounts file");

    res = 0;
err:
    free(b.obnd);
    return res;
}

int
fa_index(struct gzfh_t *fh, uint64_t optm, unsigned readlength)
{
    int len, res = -ENOMEM;
    char file[1024];
    kct_t kc = {0};
    kc.readlength = readlength;
    unsigned mode;

    C char *ext[7] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub", ".uq"};

    if (fh[0].name) {
        len = strlen(fh[0].name);
    } else {
        NB(fh[2].name != NULL);
        fh[0].name = file;

        len = strlen(fh[2].name);
        strncpy(fh[0].name, fh[2].name, ++len);
        _ACTION0(reopen(fh, ext[0], ext[5]), "");
    }
    NB(len < 256, "filename too long: %s", fh[0].name);
    if (fh[1].name == NULL)
        fh[1].name = &file[len];

    strncpy(fh[1].name, fh[0].name, len);

    if (fh[3].name == NULL)
        fh[3].name = &file[len*2];

    strncpy(fh[3].name, fh[0].name, len);

    kc.contxt_idx = buf_init_arr(kc.contxt_idx, KEYNT_BUFSZ_SHFT);
    // first check whether unique boundary is ok.
    if (fh[0].fp) {
        mode = 2;
    } else {
         _ACTION0(reopen(fh, ext[5], ext[4]), "trying to open %s instead", ext[4]);
        mode = fh[0].fp != NULL;
    }
    if (mode) {
        _ACTION(load_boundaries(fh, &kc), "loading boundary file %s", fh[0].name);
    }

    // keycount file did not exist, check twobit file.
    _ACTION0(reopen(fh, ext[mode == 2 ? 5 : 4], ext[1]),
                "%s does%s exist", ext[1], fh[0].fp ? "" : " not");
    _ACTION0(reopen(fh + 1, ext[5], ext[2]), "%s does%s exist", ext[2], fh[1].fp ? "" : " not");
    _ACTION0(reopen(fh + 3, ext[5], ext[3]), "%s does%s exist", ext[3], fh[3].fp ? "" : " not");

    if (mode && fh[0].fp && fh[1].fp && fh[3].fp) {
         _ACTION(load_seqb2(fh, &kc), "loading twobit sequence file");
         _ACTION(load_kc(fh + 3, &kc), "loading keycounts file");
    } else {
        bool found = false;
        for (int i=0; i != 4; ++i) {
            if (i == 2) continue; // i.e. fasta source!
            if (fh[i].fp) {
                if (~optm & amopt('f')) {
                    found = true;
                    EPR("%s file present, refusing to overwrite.", fh[i].name);
                }
                rclose(fh +i);
            }
        }
        if (optm & amopt('f')) {
            rclose(fh + 3+mode);
            mode = 0;
        } else if (mode) {
            EPR("%s file present, refusing to overwrite.", ext[3+mode]);
            found = true;
        }
        if (found) goto err;
        EPR("starting from scratch.");
        _ACTION(fa_read(fh, &kc), "reading fasta");
    }
    if (mode < 2) {
        _ACTION(reopen(fh, ext[1], ext[5]), "");
        _ACTION(reopen(fh + 3, ext[3], ext[6]), "");
        _ACTION(extd_uniqbnd(&kc, fh), "extending unique boundaries");
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


