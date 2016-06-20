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
    uint32_t* hkoffs = kc->hkoffs;

    for (uint32_t*k = kc->kct; k - kc->kct < kc->kct_l;) {
        while (k - kc->kct < *hkoffs) {
            char c = k!=tk?((b.prev == NO_K || k!=kc->kct + kc->contxt_idx[b.prev])?(k!=b.tgtk?' ':'t'):'p'):'k';
            EPR0("%c>%s:%u\t%u\t%c", c, kc->id + h->ido, i, *hkoffs, *k & DUP_BIT?'*':' ');
            if (*k != 0)
                print_posseq(s, *k);
            else
                EPR("(removed)");
            ++k;
        }
        // ensure we continue printing if this contig was not yet finished (for debugging)
        if (++hkoffs == kc->hkoffs + kc->hkoffs_l)
           hkoffs = &kc->kct_l;
        if (h == kc->h + kc->h_l - 1) {
            h = kc->h;
            s = kc->s;
            ++i; // increase extenion nr.
        } else {
            s += h++->len;
        }
    }
}

// if kct grows, pointers *thisk, k, b.tgtk b.prev become invalid
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
 * b.it contains ranges per contig for sequence not yet be mappable, but may become so.
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
shrink_mantra(kct_t *kc, Bnd &b, uint32_t C*C k)
{
    NB(b.it != kc->bnd->end());
    if (b.prev != NO_K) { // not at mantra start

        Mantra copy = *b.it;
        (*b.it).e = prev_pos(kc, b);
        kc->bnd->insert(b.it, copy);
    }
    (*b.it).s = b2pos_of(*k) + 2;
    //GDB:mantra1
}

static void
k_compression(kct_t *kc, Bnd &b, uint32_t *k)
{
    Hdr* h = kc->h;
    b.s = kc->s;
    b.it = kc->bnd->begin();
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
            if (b.it != kc->bnd->end() && b2pos_of(*b.tgtk) >= (*b.it).e) { // mantra end
                (*b.it).e = b2pos_of(*b.tgtk);
                ++b.it;
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
    kc->kct[kc->kct_l] = seq.p = *k;
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s, 0);
    *k ^= *k;//
    *contxt_idx = kc->kct_l++;//GDB:1
    ++b.moved;
    return contxt_idx;
}

static uint32_t *
excise(kct_t *kc, Bnd &b, uint32_t **thisk)
{
    uint32_t *contxt_idx = NULL;
    for (uint32_t *k = b.tgtk + b.moved; k < *thisk; ++k)
        contxt_idx = excise_one(kc, b, thisk, k);

    return contxt_idx;
}

static inline uint32_t *
move_uniq_one(kct_t *kc, Bnd &b, keyseq_t &seq, uint32_t *contxt_idx, C uint32_t pend)
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

        if (b2pos_of(*k) >= (*b.it).e)
            EPR("// a position is pending for %u'th, first was excised", seq.p>>1);

        *k = seq.p; // set new pos and strand, unset dupbit
        if (b.tgtk != k) {
            NB(*k);
            *b.tgtk = *k;
            *k ^= *k;//
            *contxt_idx = b.tgtk - kc->kct;//GDB:move
        }
        ++b.tgtk;

        NB(b.tgtk <= kc->kct + kc->kct_l);
    }
    seq.p = b2pos_of(seq.p);
    if (seq.p == pend)
        return NULL;

    get_next_nt_seq(b.s, seq);
    contxt_idx = get_kct(kc, seq, 1);
    seq.p += 2;
    return contxt_idx;
}

static keyseq_t
move_uniq(kct_t *kc, Bnd &b, C uint32_t pend)
{
    keyseq_t seq = {.p = after_prev(kc, b)};
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s);// already increments seq.p
    NB(*contxt_idx < kc->kct_l);
    NB(*contxt_idx != NO_K);

    while (contxt_idx)
        contxt_idx = move_uniq_one(kc, b, seq, contxt_idx, pend); //GDB

    return seq;
}

static uint32_t*
process_mantra(kct_t *kc, Bnd &b, uint32_t *k)
{
    int scope = in_scope(kc, b, k);
    uint32_t *contxt_idx = NULL;

    if (scope >= 0) { // a 2nd uniq
        shrink_mantra(kc, b, k);
        // The distance between b.tgtk and k may have grown by excised 1st kcts (beside uniq).
        excise(kc, b, &k);
        contxt_idx = excise_one(kc, b, &k, k);

        NB(contxt_idx != NULL);
        b.prev = contxt_idx - kc->contxt_idx;//GDB:moved
        EPR("previously moved were %u", b.moved);
        return k;
    }
    // prev to pend are uniques, not in scope. between uniqs are
    // from prev to p, add position if pending and reevaluate dupbit

    uint32_t tp = b2pos_of(*k);
    keyseq_t seq = move_uniq(kc, b, tp < (*b.it).e ? tp - 2 : (*b.it).e);

    seq.p = b2pos_of(seq.p);
    get_next_nt_seq(b.s, seq);
    contxt_idx = get_kct(kc, seq, 0);

    if (tp != (*b.it).e - 2) { //excise just one unique

        EPR("only one uniq isolated from mantra");
        ++b.moved; //namely last uniq.
        buf_grow_ks(kc, b, &k, NULL);
        kc->kct[kc->kct_l] = *k;
        *contxt_idx = kc->kct_l++;
        *k ^= *k;//

        Mantra copy = *b.it;
        (*b.it).s = tp + 2;
        copy.e = tp;
        kc->bnd->insert(b.it, copy);
    }
    b.prev = contxt_idx - kc->contxt_idx;//GDB:2
    return k;
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
ext_uq_iter(kct_t *kc)
{
    uint32_t *k = kc->kct;
    Bnd b = {
        .tgtk = k,
        .s = kc->s,
        .prev = NO_K,
        .moved = 0,
        .it = kc->bnd->begin()
    };
    uint32_t skctl = kc->kct_l;
    Hdr* h = kc->h;

    do {
        while (h - kc->h != (*b.it).ho) {
            // also update header
            buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
            b.s += h++->len;
        }
        b.prev = NO_K;
        uint32_t* hkoffs = kc->hkoffs + (*b.it).ho;

        NB(k - kc->kct <= *hkoffs);
        while ((k - kc->kct) < *hkoffs && b2pos_of(*k) < (*b.it).e) {

            if (IS_UQ(k)) //GDB:UQ1
                k = process_mantra(kc, b, k);

            ++k;
        }

        uint32_t end = k - kc->kct == *hkoffs ? h->end : b2pos_of(*k);
        // check whether last uniq was adjoining end
        if (end < (*b.it).s + (kc->extension << 1)) {
            excise(kc, b, &k);
            b.it = kc->bnd->erase(b.it);
        } else if (b.prev != NO_K && end < prev_pos(kc, b) + (kc->extension << 1)) {
            excise(kc, b, &k);
            (*b.it).e = prev_pos(kc, b);
            ++b.it;
        } else {

            move_uniq(kc, b, end);

            ++b.it;
        }
        *hkoffs = b.tgtk - kc->kct;
        //GDB:next mantra

    } while (b.it != kc->bnd->end());

    k_compression(kc, b, k);//GDBk
    NB(skctl == kc->kct_l); //GDBk
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res = -ENOMEM;
    for (kc->extension = 0; kc->extension != kc->readlength - KEY_WIDTH + 1; ++kc->extension) {
        kc->iter = 0;
        do { // until no no more new uniques
            kc->uqct = 0;
            ext_uq_iter(kc);
            EPR("observed %u excised in iteration %u, extension %u\n",
                kc->uqct, ++kc->iter, kc->extension);
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


