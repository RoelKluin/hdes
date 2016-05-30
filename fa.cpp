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
        while (k - kc->kct <= *hkoffs) {
            char c = k!=tk?(k!=b.prev?(k!=b.tgtk?' ':'s'):'p'):'k';
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

/*
 * b.it contains ranges per contig for sequence not yet be mappable, but may become so.
 * Initially there are one or more dependent on the presence and number of N-stretches;
 * initialized in key_init.cpp, But as adjacent uniques cover regions, the remaining `mantra'
 * shrinks (here).
 *
 * mantras start at the chr start, the position after a N-stretch or after a region scoped
 * by uniques, the first non-unique after. The mantra ends at chr end, stretch or first non
 * unique before a region scoped by uniques.
 */
static void
shrink_mantra(kct_t *kc, Bnd &b, uint32_t C*C k)
{
    NB(b.it != kc->bnd->end());
    if (b.prev) { // not at mantra start
        C uint32_t ke = (*b.it).ke;
        (*b.it).ke = b.prev - kc->kct;

        kc->bnd->insert(b.it, *b.it);     // copy of current
        // mantra became two smaller ranges.
        (*b.it).ke = ke;
    }
    (*b.it).s = b2pos_of(*k) + 2;
    //GDB:mantra1
}

// if kct grows, pointers *thisk, k, b.tgtk b.prev become invalid
static void
buf_grow_ks(kct_t *kc, Bnd &b, uint32_t **k1, uint32_t **k2)
{
    if ((kc->kct_l + 1) >= (1ul << kc->kct_m)) {
        uint32_t ok1, ok2, ok3, ok4;

        ok1 = b.tgtk - kc->kct;
        // adapting b.prev is not needed: is always set after function.
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

static uint32_t *
excise(kct_t *kc, Bnd &b, uint32_t **thisk)
{
    uint32_t *contxt_idx = NULL;
    keyseq_t seq = {0};
    for (uint32_t *k = b.tgtk + b.moved; k <= *thisk; ++k) {
        NB(k < kc->kct + kc->kct_l);
        buf_grow_ks(kc, b, thisk, &k);
        NB(k != kc->kct + kc->kct_l);
        kc->kct[kc->kct_l] = seq.p = *k;
        contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s, 0);
        *k ^= *k;//
        *contxt_idx = kc->kct_l++;//GDB:1
        ++b.moved;
    }
    return contxt_idx;
}

static void
process_mantra(kct_t *kc, Bnd &b, uint32_t **thisk)
{

    uint32_t *contxt_idx = NULL;

    if (in_scope(kc, b, *thisk)) { // a 2nd uniq
        shrink_mantra(kc, b, *thisk);
        // The distance between b.tgtk and k may have grown by excised 1st kcts (beside uniq).
        contxt_idx = excise(kc, b, thisk);
        NB(contxt_idx != NULL);
        b.prev = kc->kct + *contxt_idx;
        EPR("previously moved were %u", b.moved);
        return;
    }
    // prev to pend are uniques, not in scope. between uniqs are
    // from prev to p, add position if pending and reevaluate dupbit
    keyseq_t seq = { .p = after_prev(b) };
    contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s); // already increments seq.p
    NB(*contxt_idx != NO_K);

    C uint32_t pend = before_this(kc, b, *thisk);
    for (;;) {
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

            if (k - kc->kct > (*b.it).ke)
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
            break;

        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq, 1);
        seq.p += 2;
    }
    if (*thisk - kc->kct != (*b.it).ke) { //excise just one unique

        EPR("only one uniq isolated from mantra");
        ++b.moved; //namely last uniq.
        seq.p = b2pos_of(seq.p);
        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq, 0);
        buf_grow_ks(kc, b, thisk, NULL);
        uint32_t *k = kc->kct + kc->kct_l;

        NB(k != *thisk);
        *k = **thisk;
        b.prev = k;
        shrink_mantra(kc, b, *thisk);
        **thisk ^= **thisk;//
        *contxt_idx = kc->kct_l++;//GDB:2

    } else {
        b.prev = *thisk;
    }
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
static int
ext_uq_iter(kct_t *kc)
{
    uint32_t *k = kc->kct; // location after uniques, from last time
    Bnd b = {
        .tgtk = k,
        .s = kc->s,
        .prev = NULL,
        .moved = 0,
        .it = kc->bnd->begin()
    };
    uint32_t skctl = kc->kct_l;
    Hdr* h = kc->h;
    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)

    while (1) {

        if (IS_UQ(k)) // they are handled during excision (or postponed) - process_mantra()
            process_mantra(kc, b, &k);

        if (k - kc->kct == (*b.it).ke) {
            if (k == hdr_end_k(kc, h)) {
                // check whether last uniq was adjoining end
                if (b.prev) {
                    // in scope ?
                    if (h->end <= (kc->extension << 1) + b2pos_of(*b.prev)) {
                        EPR("in scope...");
                        (*b.it).ke = b.prev - kc->kct;
                        excise(kc, b, &k);
                    } else {
                        EPR("// XXX single uniq excision (or was it already done?)");
                    }
                } else {
                    b.it = kc->bnd->erase(b.it); // nothing left
                }

                buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l - 1);
                // also update new end for header
                kc->hkoffs[h - kc->h] = b.tgtk - kc->kct - 1;

                if (b.it == kc->bnd->end())
                    break;

                b.s += h->len;
                ++h;
            }
            b.prev = NULL;

            // next mantra
            (*b.it).ke = b.tgtk - kc->kct - 1;
            if (++b.it == kc->bnd->end())
                break;
        }
        ++k;
    }

    kc->uqct = k - b.tgtk;
    h = kc->h;
    b.s = kc->s;
    uint32_t* hkoffs = kc->hkoffs + kc->hkoffs_l - kc->h_l;
    // put uniqs and 1st excised back in array (recompression)
    for (k = kc->kct + skctl; k - kc->kct < kc->kct_l; ++h, ++hkoffs) {

        while (k - kc->kct <= *hkoffs) {
            if (*k != 0) {
                keyseq_t seq = {.p = *k };
                kc->contxt_idx[build_ndx_kct(kc, seq, b.s, 0)] = b.tgtk - kc->kct;
                *b.tgtk++ = *k;
                *k ^= *k;//
            }
            ++k;
        }
        NB(hkoffs < kc->hkoffs + kc->hkoffs_l);
        *hkoffs = b.tgtk - kc->kct - 1;
        b.s += h->len;
    }
    kc->kct_l = skctl;
    return 0;
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res = -ENOMEM;
    for (kc->extension = 1; kc->extension != kc->readlength - KEY_WIDTH + 1; ++kc->extension) {
        kc->iter = 0;
        do { // until no no more new uniques
            kc->uqct = 0;
            _EVAL(ext_uq_iter(kc));
            EPR("observed %u uniques in iteration %u, extension %u\n",
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


