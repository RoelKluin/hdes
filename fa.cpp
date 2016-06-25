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
#include <limits.h> //INT_MAX
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <set>
#include <algorithm> // swap
#include <setjmp.h>
#include <cmocka.h> // TODO: unit testing
//#include <stack>
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"
#include "b6.h"

void
free_kc(kct_t *kc)
{
    delete kc->bnd;
    buf_free(kc->id);
    buf_free(kc->s);
    buf_free(kc->kct);
    buf_free(kc->h);
    buf_free(kc->hkoffs);
    buf_free(kc->contxt_idx);
}

static void
print_kct(kct_t *kc, Bnd &b, uint32_t* tk)
{
    unsigned i = 0;
    EPR("");
    Hdr* h = kc->h;
    uint8_t *s = kc->s;
    std::list<Mantra>::iterator it = kc->bnd->begin();
    uint32_t* hkoffs = kc->hkoffs;

    for (uint32_t*k = kc->kct; k - kc->kct < kc->kct_l;) {
        if (it != kc->bnd->end()) {
            while (h - kc->h != (*it).ho)
                s += h++->len;
            hkoffs = kc->hkoffs + (*it).ho;
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
        if (it != kc->bnd->end())
            ++it;
    }
}

// if kct grows, pointers *thisk, k, b.tgtk become invalid
static void
buf_grow_ks(kct_t *kc, Bnd &b, uint32_t **k1, uint32_t **k2)
{
    if ((kc->kct_l + 1) >= (1ul << kc->kct_m)) {
        uint32_t ok1, ok3, ok4;

        ok1 = b.tgtk - kc->kct;
        ok3 = *k1 - kc->kct;
        ok4 = k2 ? *k2 - kc->kct : ~0u;
        uint32_t *t = (uint32_t *)realloc(kc->kct, sizeof(uint32_t) << ++kc->kct_m);
        if_ever (t == NULL)
            raise(SIGTRAP);
        kc->kct = t;
        b.tgtk = t + ok1;
        *k1 = t + ok3;
        if (ok4 != ~0u) *k2 = t + ok4;
    }
}

/*
 * kc->bnd contains ranges per contig for sequence not yet be mappable, but may become so.
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
k_compression(kct_t *kc, Bnd &b, uint32_t *k)
{
    Hdr* h = kc->h;
    b.s = kc->s;
    std::list<Mantra>::iterator it = kc->bnd->begin();
    uint32_t *hkoffs = kc->hkoffs + kc->h_l;
    // put uniqs and 1st excised back in array (recompression)
    while(1) {
        for (;k - kc->kct < *hkoffs; ++k) {

            if (*k == 0)
                continue;

            if (k != b.tgtk) {
                keyseq_t seq = {.p = *k };
                kc->contxt_idx[build_ndx_kct(kc, seq, b.s, 0)] = b.tgtk - kc->kct;
                *b.tgtk = *k;
                *k ^= *k;//
            }
            if (it != kc->bnd->end() && b2pos_of(*b.tgtk) >= (*it).e) { // mantra end
                (*it).e = b2pos_of(*b.tgtk);
                ++it;
            }
            ++b.tgtk;
        }
        *hkoffs = b.tgtk - kc->kct;
        b.s += h->len;
        if(++h == kc->h + kc->h_l) {
            h = kc->h;
            b.s = kc->s;
        }
        if (hkoffs == &kc->kct_l)
            break;

        if (++hkoffs == kc->hkoffs + kc->hkoffs_l)
            hkoffs = &kc->kct_l;
    }
    kc->kct_l = b.tgtk - kc->kct;
    buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
}

static inline uint32_t *
excise_one(kct_t *kc, Bnd &b, uint32_t **thisk, uint32_t *k)
{
    keyseq_t seq = {0};
    NB(k < kc->kct + kc->kct_l);
    buf_grow_ks(kc, b, thisk, &k);
    NB(k != kc->kct + kc->kct_l);
    //*k &= ~DUP_BIT; // or keys that were moved after extension still have their dup bit set.
    kc->kct[kc->kct_l] = seq.p = *k;
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s, 0);
    *k ^= *k;//
    *contxt_idx = kc->kct_l++;
    ++b.moved;
    return contxt_idx;//P;
}

static uint32_t *
excise(kct_t *kc, Bnd &b, uint32_t **thisk)
{
    uint32_t *contxt_idx = NULL;
    for (uint32_t *k = b.tgtk + b.moved; k < *thisk; ++k)
        contxt_idx = excise_one(kc, b, thisk, k);

    return contxt_idx;
}

static inline void
move_uniq_one(kct_t *kc, Bnd &b, keyseq_t &seq, uint32_t *contxt_idx)
{
    NB(*contxt_idx != NO_K);
    NB(*contxt_idx < kc->kct_l);
    uint32_t *k = kc->kct + *contxt_idx;
    print_seq(&seq);

    if (k < b.tgtk) {
        // second occurance

        if (~*k & DUP_BIT) {
            // no dup after all

            *k |= DUP_BIT;
            --kc->uqct;
        }
    } else {

        // 1st occurance;
        ++kc->uqct;    // unique or decremented later.

        *k = seq.p; // set new pos and strand, unset dupbit XXX:Invalid write of size 4
        if (b.tgtk != k) {
            NB(*k);
            *b.tgtk = *k;
            *k ^= *k;//
            *contxt_idx = b.tgtk - kc->kct;//K; moved up
        }
        ++b.tgtk;

        NB(b.tgtk <= kc->kct + kc->kct_l);
    }
}

static keyseq_t
move_uniq(kct_t *kc, Bnd &b, C uint32_t start, C uint32_t pend)
{
    keyseq_t seq = {.p = start};
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s);// already increments seq.p
    seq.p = b2pos_of(seq.p);
    NB(*contxt_idx < kc->kct_l);
    NB(*contxt_idx != NO_K);
    NB(seq.p <= pend);

    while (seq.p != pend) {
        move_uniq_one(kc, b, seq, contxt_idx);
        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq, 1);
        seq.p = b2pos_of(seq.p) + 2;
    }
    move_uniq_one(kc, b, seq, contxt_idx);

    return seq;
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
ext_uq_iter(kct_t *kc, uint32_t ext)
{
    uint32_t *k = kc->kct;
    Bnd b = {
        .tgtk = k,
        .s = kc->s,
        .moved = 0,
    };
    std::list<Mantra>::iterator it = kc->bnd->begin();
    Hdr* h = kc->h;//B; initial state
    uint32_t skctl = kc->kct_l;//B; initial state

    do {
        NB(h - kc->h <= (*it).ho);
        while (h - kc->h != (*it).ho) {
            //~ also update header
            buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
            b.s += h++->len;
        }
        uint32_t* hkoffs = kc->hkoffs + (*it).ho;

        //k - kc->kct <= *hkoffs: may be untrue after key excision.
        uint32_t end = (*it).e;

        while ((k - kc->kct) < *hkoffs && b2pos_of(*k) < end) {

            if (IS_UQ(k)) { //~ uniq
                end = b2pos_of(*k);
                break;
            }
            ++k;
        }

        if ((*it).s + ext >= end) {
            if (end != (*it).e) {
                (*it).s = end + 2;
                ++k;
            } else {
                it = kc->bnd->erase(it);
                *hkoffs = b.tgtk - kc->kct;
            }
            excise(kc, b, &k);
        } else {
            move_uniq(kc, b, (*it).s, end - 2);
            if (end != (*it).e) {
                //if (end + 2 != (*it).e) {
                Mantra copy = *it;
                copy.e = end;
                (*it).s = end + 2;
                kc->bnd->insert(it, copy);
                ++k;
            } else {
                ++it;
                *hkoffs = b.tgtk - kc->kct;
            }
        }

    } while (it != kc->bnd->end());//B; next region

    k_compression(kc, b, k);
    NB(skctl == kc->kct_l);//B; final state
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res = -ENOMEM;
    uint32_t end = (kc->readlength - KEY_WIDTH + 1) << 1;
    for (uint32_t ext = 0; ext != end; ext += 2) {
        kc->iter = 0;
        do { // until no no more new uniques
            kc->uqct = 0;
            ext_uq_iter(kc, ext);
            EPR("observed %u excised in iteration %u, extension %u\n",
                kc->uqct, ++kc->iter, ext);
        } while (kc->uqct > 0);
    }
    _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
    _ACTION(save_kc(fhout + 3, kc), "writing unique keycounts file");
err:
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


