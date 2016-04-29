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
//#include <stack>
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"
#include "b6.h"

void
free_kc(kct_t *kc)
{
    for (Hdr *h = kc->h; h != kc->h + kc->h_l; ++h) {
        free(h->part);
        delete h->bnd;
    }
    _buf_free(kc->id);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->h);
    _buf_free(kc->hk);
    _buf_free(kc->ndxkct);
}

static void
swap_kct(seq_t *kcnxk,  uint64_t *k1,  uint64_t *k2, seq_t *ndxkct2, uint8_t C*C s)
{
    if (k1 == k2)
        return;
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {.p = b2pos_of(*k1)};
    seq_t *ndxkct1 = kcnxk + build_ndx_kct(seq, s);
    std::swap(*ndxkct1, *ndxkct2); // ndxkct points to reordered kcts
    std::swap(*k1, *k2);           // and swap info at kcts
}

/*
 * b.it contains ranges per contig for sequence not yet be mappable, but may become so.
 * Initially there are one or more dependent on the presence and number of N-stretches;
 * initialized in key_init.cpp, But as adjacent uniques cover regions, the remaining `mantra'
 * shrinks (here).
 */
static void
unique_covered(Bnd &b, uint64_t C*C thisk, C pos_t prev, C pos_t p)
{
    if (b.prev) {                          // if not at boundary start
        C pos_t end = (*b.it).e;           // store original end
        (*b.it).e = prev;                  // shift end
        if (thisk) {                       // if not at boundary end either
            b.h->bnd->insert(b.it, *b.it); // insert a copy of the current
            (*b.it).s = p;                 // mantra became divided in two smaller ranges.
            (*b.it).e = end;               // reinstate original end
        }                                  // (otherwise the last boundary just got smaller)
    } else if (thisk) {     // if at chr start, but not at chr end
        (*b.it).s = p;      // shift start
    } else {                              // entire region became mapable
        b.it = b.h->bnd->erase(b.it);     // this returns the next element
        b.prev = NULL;
    }
}

/*
 * If two unique keys occur within scope, all reads (i.e. KEY_WIDTH + current extension)
 * within these are mapable due to the two. This site no longer needs to be considered
 * for the inbetween - these keys may henceforth become secondary uniques elsewhere.
 *
 * Of unique keys the stored position is guaranteed correct. (of non-uniques there is only
 * the garuantee that the stored position results in the key that corresponds).
 *
 */
static inline void
extd_uq_by_p(kct_t *kc, pos_t p, pos_t pend, uint8_t C*C s, uint64_t C*C sk)
{
    NB(pend < (kc->s_l << 2) && p < pend);

    if (p + kc->ext + 1 < pend || p + 1 == pend) // do we have keys between unique?
        return;                                  // no: then nothing to do here.

    keyseq_t seq = {.p = p + 1};
    seq_t *ndxkct = kc->ndxkct + build_ndx_kct(seq, s);

    for(;;) {
        uint64_t *k = kc->kct + *ndxkct;

        // stored first occurance matches current position and contig
        if (b2pos_of(kc, k) == seq.p && k >= sk) {

            NB(*k & DUP_BIT);    // the first occurance of the key is excised.
            ++kc->reeval;        // leave dup bit set to mark a position is pending
        }
        if (++seq.p == pend)
            break;

        get_next_nt_seq(s, seq);
        ndxkct = get_kct(kc, seq);
    }
}

static void
reached_boundary(kct_t *kc, Bnd &b)
{
    C pos_t prev = _prev_or_bnd_start(b);
    if ((prev + kc->ext > (*b.it).e) && ((*b.it).e > prev)) {
        if (b.prev) {
            ++b.it;
        } else {
            b.it = b.h->bnd->erase(b.it);
        }
        extd_uq_by_p(kc, (*b.it).s, prev, kc->s + b.h->s_s, b.sk);
    }
}



static void
handle_range(kct_t *kc, Bnd &b, uint64_t C*C thisk)
{
    uint8_t C*C s = kc->s + b.h->s_s;
    pos_t prev = _prev_or_bnd_start(b);
    C pos_t pend = thisk ? b2pos_of(*thisk) : (*b.it).e;
    NB(prev < pend);

    if (pend - prev <= kc->ext + 1) { // a 2nd uniq in scope

        unique_covered(b, thisk, prev, pend);
        extd_uq_by_p(kc, prev, pend, s, b.sk);
        return;

    }

    // first uniq in scope. from prev to p, add position if pending and reevaluate dupbit
    keyseq_t seq = { .p = prev + 1};
    seq_t *ndxkct = kc->ndxkct + build_ndx_kct(seq, s);

    for(;;) {
        uint64_t *k = kc->kct + *ndxkct;

        if (~*k & DUP_BIT) {

            if (k < b.sk) {
                *k |= DUP_BIT;      // second occurance, no dup after all
                --kc->reeval;
            } else {

                EPR("uniq");
                print_seq(&seq);
                *k |= DUP_BIT; // FIXME don't want this.
                //swap_kct(kc->ndxkct, b.sk++, k, ndxkct, s);
            }

        } else if (k >= b.sk) { // 1st occurance, place it now.

            if (b2pos_of(kc, k) == seq.p) {

                ++kc->reeval;   // may yet be unique or we decrement this later.

            } else {
                // a position is pending,  see extd_uq_by_p()

                NB(b2pos_of(kc, k) < seq.p || b.hk == kc->hk || k < kc->kct + (b.hk-1)->koffs);

                *k &= INFIOR_MASK; // unset strand bit and pos, dupbit (highest) is kept.
                *k |= (uint64_t)!seq.t << ORIENT_SHFT | seq.p; // set new pos and strand
            }

            *k &= ~DUP_BIT;
            swap_kct(kc->ndxkct, b.sk++, k, ndxkct, s);
        }
        if (++seq.p == pend)
            break;

        get_next_nt_seq(s, seq);
        ndxkct = get_kct(kc, seq);
    }
}

static void
update_header(kct_t *kc, uint64_t *k, Bnd &b)
{
    while (k >= kc->kct + b.hk->koffs) {
        handle_range(kc, b, NULL);
        reached_boundary(kc, b);
        NB(b.hk < kc->hk + kc->hk_l);
        b.h = kc->h + (++b.hk)->hoffs;
        b.it = b.h->bnd->begin();
        b.prev = NULL;
    }
}

static void
update_boundary(kct_t *kc, uint64_t *k, Bnd &b)
{
    while (b2pos_of(*k) >= (*b.it).e) {
        handle_range(kc, b, NULL);
        reached_boundary(kc, b);
        ++b.it;
        NB(b.it != b.h->bnd->end());
        b.prev = NULL;
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
static void
ext_uq_iter(kct_t *kc)
{
    uint64_t *k = kc->kct;
    uint64_t *kend = k + kc->kct_l - 1 - kc->last_uqct; // location after uniques, from last time
    Bnd b = {
        .sk = k,
        .hk = kc->hk,
        .h = kc->h + kc->hk->hoffs,
        .prev = NULL,
        .it = kc->h->bnd->begin()
    };

    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    for (; k != kend; ++k) {

        if (IS_DUP(k)) // they are handled during excision (or postponed) - handle_range()
            continue;

        update_header(kc, k, b);
        update_boundary(kc, k, b);
        handle_range(kc, b, k);
        b.prev = k;
    }
    handle_range(kc, b, NULL);
    reached_boundary(kc, b);

    //NB(0, "(NOT!;) end of loop reached! (TODO: swapping..)");
    // FIXME: iteration over swapped instead of uniques
    //

    // TODO: some uniques no longer occur: all are excised.
    kc->last_uqct = kc->uqct;

    EPQ(dbg > 0, "observed %u uniques in iteration %u, extension %u, %u to be reevaluated\n",
            kc->uqct, ++kc->iter, kc->ext, kc->reeval);
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res;
    kc->uqct = kc->reeval = 0;

    for (kc->ext = res = 1; res > 0 && kc->ext != kc->readlength - KEY_WIDTH + 1;++kc->ext) {
        kc->iter = 0;
        do { // until no no more new uniques
            ext_uq_iter(kc);
        } while (kc->reeval > 0);
    }
    if (res == 0) {
        _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
        _ACTION(save_kc(fhout + 3, kc), "writing unique keycounts file");
    }
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

    kc.ndxkct = _buf_init_arr_err(kc.ndxkct, KEYNT_BUFSZ_SHFT, return -ENOMEM);
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


