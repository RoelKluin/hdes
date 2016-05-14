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
    _buf_free(kc->id);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->h);
    _buf_free(kc->hk);
    _buf_free(kc->ext);
    _buf_free(kc->contxt_idx);
}

static void
swap_kct(uint32_t *kcnxk,  uint32_t *k1,  uint32_t *k2, uint32_t *contxt_idx2, uint8_t C*C s)
{
    if (k1 == k2)
        return;
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {.p = b2pos_of(*k1)};
    uint32_t *contxt_idx1 = kcnxk + build_ndx_kct(seq, s);
    std::swap(*contxt_idx1, *contxt_idx2); // contxt_idx points to reordered kcts
    std::swap(*k1, *k2);           // and swap info at kcts
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
shrink_mantra(kct_t *kc, Bnd &b, uint32_t C*C thisk, C uint32_t prev, C uint32_t p)
{
    if (b.prev) {                          // not at mantra start
        C uint32_t end = (*b.it).e;
        (*b.it).e = prev;
        if (thisk) {                       // not at mantra end either
            kc->bnd->insert(b.it, *b.it);  // copy of current
            (*b.it).s = p;                 // mantra became two smaller ranges.
            (*b.it).e = end;
        } // or end got shifted
    } else {
        // shift start or if at end, force erase (entirely mapable).
        (*b.it).s = thisk ? p : (*b.it).e; // or entire region became mapable, force erase.
    }
    NB((*b.it).s <= (*b.it).e);
    if ((*b.it).s == (*b.it).e) {
        b.it = kc->bnd->erase(b.it);   // nothing left
        b.prev = NULL;
    }
}

/*This function swaps d elements starting at index fi
  with d elements starting at index si */
void swap(uint32_t *a, uint32_t *b, uint32_t C*C c)
{
    while (a != c) {
        uint32_t temp = *a;
        *a++ = *b;
        *b++ = temp;
    }
}

/*Function to left rotate s[] of siz n by d
 * Time complexity: O(n)
 */
void leftRotate(uint32_t *a, unsigned C d, unsigned n)
{
    if(d == 0 || d == n)
        return;
    uint32_t *C b = a + d;
    uint32_t *c = a + n - d;
    while (c != b) {

        if(c >= b) {
            swap(a, c, b);
            c -= b - a;
        } else {
            swap(a, b, c);
            n = c - a;
            a = c;
            c += n;
        }
    }
    swap(a, b, b);
}


static void
process_mantra(kct_t *kc, Bnd &b, uint32_t *C thisk)
{
    C uint32_t prev = _prev_or_bnd_start(b);
    C uint32_t pend = thisk ? b2pos_of(*thisk) - 1 : (*b.it).e;

    if (in_scope(kc, prev, pend)) { // a 2nd uniq
EPR("excision");
        shrink_mantra(kc, b, thisk, prev, pend);
        return;
    }
    // prev to pend are uniques, not in scope. between uniqs are
    // from prev to p, add position if pending and reevaluate dupbit
    keyseq_t seq = { .p = prev + 1};
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(seq, b.s); // already increments seq.p
    NB(*contxt_idx != NO_KCT);
EPR("kept %u-%u", prev, pend);

    // rotation must be undone before we add to b.sk.
    // this should also fix the ndx => wrong kct_i that occurs after swaps.
    leftRotate(thisk, b.rot, b.sk - thisk);
    b.rot = 0;
    for (;;) {
        uint32_t *k = kc->kct + *contxt_idx;

        if (k < b.sk) { // XXX requires multi-contig uniqs and excised in b.sk .. k
            // second occurance

            if (~*k & DUP_BIT) {
                //no dup afterr all

                *k |= DUP_BIT;
                --kc->uqct;
            }
        } else {
            // 1st occurance
            ++kc->uqct;    // unique or decremented later.
            if (b2pos_of(kc, k) == seq.p) {

                *k &= ~DUP_BIT;

            } else {
                EPR("// a position is pending for %u'th, first was excised", seq.p);

                *k = seq.p << 1 | (seq.t != 0); // set new pos and strand, unset dupbit
            }
            ++b.rot;
            swap_kct(kc->contxt_idx, b.sk++, k, contxt_idx, b.s);
            //gdb:swap
        }
        if (seq.p == pend)
            break;

        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq);
        ++seq.p;
    }
}

static void
reached_boundary(kct_t *kc, Bnd &b)
{
    C uint32_t prev = _prev_or_bnd_start(b);

    // if at start and entirely mapable: remove mantra.
    if (b.prev == NULL && in_scope(kc, prev, (*b.it).e)) {
EPR("mantra removed");
        b.it = kc->bnd->erase(b.it);
    } else {
EPR("next bnd");
        (*b.it).ke = b.sk - kc->kct;
        ++b.it;
    }
}

static inline void
skip_mantra(kct_t *kc, Bnd &b, uint32_t *k)
{
    process_mantra(kc, b, NULL);
    reached_boundary(kc, b);
    b.prev = NULL;
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
        .sk = k,
        .s = kc->s,
        .prev = NULL,
        .rot = 0,
        .it = kc->bnd->begin()
    };

    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    for (HK *hk = kc->hk; hk != kc->hk + kc->hk_l; ++hk) {

        if (k < kc->kct + hk->koffs) {

            while (k < kc->kct + hk->koffs) {

                while (k - kc->kct < (*b.it).ke) {
EPR("%u", k - kc->kct);

                    if (IS_UQ(k)) { // they are handled during excision (or postponed) - process_mantra()
EPR("uniq");
print_posseq(b.s, b2pos_of(*k));
                        process_mantra(kc, b, k); //
                        b.prev = k;
                    }
                    ++k;
                }
                if (b.it == kc->bnd->end())
                    goto out;

                skip_mantra(kc, b, k);
            }
        } else {
            if (b.it == kc->bnd->end())
                goto out;
            skip_mantra(kc, b, k);
        }
EPR("next hdr %u, %u", hk->koffs, b.sk - kc->kct);
        NB(hk->koffs <= kc->kct_l);
        NB(hk->koffs >= b.sk - kc->kct);
        b.s += (hk->len >> 2) + !!(hk->len & 3);

        //XXX: this is cumulative for contigs for this extension and iteration.
        uint32_t uq_and_1stexcised = hk->koffs - (b.sk - kc->kct);

        _buf_grow_add_err(kc->ext, 1ul, 0, uq_and_1stexcised, return -ENOMEM);
        EPQ(uq_and_1stexcised, "added %u uniq", uq_and_1stexcised);
        kc->uqct += uq_and_1stexcised;
    }
out:
    NB(kc->uqct == k - b.sk, "%u != %u",  kc->uqct, k - b.sk);
    kc->uqct = k - b.sk;
    kc->last_uqct = kc->uqct;
    return 0;
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res = -ENOMEM;
    kc->ext = _buf_init_err(kc->ext, 1, goto err);
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

    kc.contxt_idx = _buf_init_arr_err(kc.contxt_idx, KEYNT_BUFSZ_SHFT, return -ENOMEM);
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


