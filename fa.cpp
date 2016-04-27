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

static void
start_dbg(kct_t *C kc, Hdr *C h, seq_t C dna)
{
    char C *hdr = kc->id + h->part[0];

    // TODO: if not start of contig or just after a N-stretch, the dna key should be
    // just after an unique key. we could check this.

    if (dbg == 3 || dbg == 5)
        dbg = strncmp(hdr, dbgchr, strlen(hdr)) ? 3 : 5;
    dbg = 7;

    if (dbg > 3) {
        EPR0("----[\t%s%s:(" Pfmt "+)" Pfmt ".." Pfmt "\text:%u\t]----\tdna:",
                strlen(hdr) < 8 ? "\t": "", hdr, (*kc->bdit).corr, (*kc->bdit).s,
                (*kc->bdit).e, kc->ext);
        print_dna(dna);
        show_mantras(kc, h);
    }
}

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

// reset and remove non-uniqs, move uqs to start of kc->kct array
static inline void
extd_uq_by_p(kct_t *kc, pos_t p, pos_t pend, uint8_t C*C s, C uint32_t koffs, std::set<uint64_t> &pk)
{
    seq_t *ndxkct;
    keyseq_t seq = {0};
    NB(pend < (kc->s_l << 2) && p < pend);

    if (p + kc->ext + 1 >= pend && p + 1 != pend) {
        seq.p = p;
        ndxkct = kc->ndxkct + build_ndx_kct(seq, s);
        uint64_t*k = kc->kct + *ndxkct;
        NB(kc->kct[*ndxkct] & DUP_BIT);

        for (;seq.p != pend; ++seq.p) {
            get_next_nt_seq(s, seq);
            ndxkct = get_kct(kc, seq);
            k = kc->kct + *ndxkct;

            // in next iteration: reeval whether this is still a dup.
            if (*k & DUP_BIT) {
                if (k - kc->kct < koffs || b2pos_of(kc, k) < seq.p) {
                    ++kc->reeval;
                } else if (b2pos_of(kc, k) == seq.p) {
                    // all remaining were excised this iteration
                    --kc->reeval;
                    // this position should be recognisable since it is the only with
                    // k before the adjacent uniq's k.
                    *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
                    *k |= ((uint64_t)!seq.t << ORIENT_SHFT) | seq.p;
                    pk.insert(*k);
                }
            } else { // may be unique - next iter.
                ++kc->reeval;
            }
        }
        // need to determine kcts for inbetween. Same can occur twice, use end instead.
    }
}

static void
reached_boundary(kct_t *kc, Hdr *h, uint64_t C*C last, C uint32_t koffs, std::set<uint64_t> &pk)
{
    C pos_t p = last ? b2pos_of(*last) : (*kc->bdit).e;
    show_mantras(kc, h);
    if ((p > (*kc->bdit).s) && (p - kc->ext <= (*kc->bdit).s)) {
        if (last) {
            --kc->bdit;
        } else {
            kc->bdit = h->bnd->erase(kc->bdit);
        }
        extd_uq_by_p(kc, (*kc->bdit).s, p, kc->s + h->s_s, koffs, pk);
    }
}

static void
swap_kct(kct_t *kc,  uint64_t *k1,  uint64_t *k2, seq_t *ndxkct2, uint8_t C*C s)
{
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {.p = b2pos_of(kc, k1)};
    seq_t *ndxkct1 = kc->ndxkct + build_ndx_kct(seq, s);
    std::swap(*ndxkct1, *ndxkct2); // ndxkct points to reordered kcts
    std::swap(*k1, *k2);           // and swap info at kcts
}

// XXX XXX hk need insertions.
static void
place_uniques(kct_t *kc, uint64_t *sk, uint64_t* &kend, uint8_t C*C s, std::set<uint64_t> &pk)
{
    HK *hk = kc->hk;
    NB(0, "FIXME: hk need insertions. koffs need to be updated");

    for (std::set<uint64_t>::iterator it = pk.begin(); it != pk.end(); ++it) {
        uint64_t kv = *it;
        if (~kv & DUP_BIT) {
            ++kc->uqct;
            keyseq_t seq = {.p = b2pos_of(kv)};
            seq_t *ndxkct = kc->ndxkct + build_ndx_kct(seq, s);
            // If smaller than sk, an ambiguous key was apparently unique. These are handled in
            // next iteration
            NB(kc->kct + *ndxkct >= sk);

            swap_kct(kc, kc->kct + *ndxkct, kend--, ndxkct, s);
        }
    }
}

/*
 * kc->bdit contains ranges per contig for sequence cannot yet be mapped, but may become so.
 * Initially there are one or more dependent on the presence and number of N-stretches.
 * initialized in key_init.cpp, But as adjacent uniques cover regions, the remaining `mantra'
 * shrinks (here).
 */
static void
uq_bounded(kct_t *kc, Hdr *h, uint64_t C*C prevk, uint64_t C*C thisk, C pos_t prev, C pos_t p)
{
    if (prevk) { NB(prev < p);
        C pos_t b2end = (*kc->bdit).e;
        (*kc->bdit).e = prev; // shift end
        if (thisk) {
            h->bnd->insert(kc->bdit, *kc->bdit); // insert a copy of the current
            (*kc->bdit).e = b2end; // restore orig end
        }
    } else {
        (*kc->bdit).s = p; // shift start
    }
}


static void
handle_range(kct_t *kc, Hdr *h, uint64_t C*C prevk, uint64_t C*C thisk,
        C uint32_t koffs, uint64_t** sk, std::set<uint64_t> &pk)
{
    seq_t *ndxkct;
    keyseq_t seq = {0};
    uint8_t C*C s = kc->s + h->s_s;
    C pos_t b2start = (*kc->bdit).s;
    C pos_t b2end = (*kc->bdit).e;
    pos_t prev = prevk ? b2pos_of(*prevk) : b2start + KEY_WIDTH - 1;
    C pos_t p = thisk ? b2pos_of(*thisk) : b2end; NB(prev < p);

    if (p - prev <= kc->ext + 1) { // a 2nd uq
        uq_bounded(kc, h, prevk, thisk, prev, p);

        // what if all remaining occurances happen between two boundaries?
        if (prev - p > 1)
            extd_uq_by_p(kc, prev, p, s, koffs, pk);
        seq.p = prev = p;
        ndxkct = kc->ndxkct + build_ndx_kct(seq, s);

    } else { // first occurance.

        // 1: from prev till now, add position and dup reevaluation
        seq.p = prev - KEY_WIDTH + 1;
        build_key(s, seq, prev + 1);

        do {
            ndxkct = get_kct(kc, seq);
            uint64_t *k = kc->kct + *ndxkct;

            if (*k & DUP_BIT) {
                if (k - kc->kct < koffs || b2pos_of(kc, k) < seq.p) {
                    *k &= ~DUP_BIT; NB(*sk <= k); // unset for 1st occurance
//                    swap_kct(kc, *sk, k, ndxkct, s); //XXX XXX: valid ???
                    ++kc->reeval;
//                    k = (*sk)++;
                } else if (b2pos_of(*k) == seq.p) {
                    // all but last positions were excised, unique after all!
                    uq_bounded(kc, h, prevk, thisk, prev, seq.p);
                    pk.insert(*k);
                    unsigned at = seq.p;
                    if (at < prev + kc->ext + 1) {
                        if (at > prev + 1)
                            extd_uq_by_p(kc, prev, at, s, koffs, pk);
                        prev = at;
                        if (p - prev < kc->ext + 1) {
                            //rollback entirely
                            uq_bounded(kc, h, prevk, NULL, prev, p);
                            seq.p = p;
                            ndxkct = kc->ndxkct + build_ndx_kct(seq, s);
                            if (p - prev > 1)
                                extd_uq_by_p(kc, prev, p, s, koffs, pk);
                            prev = p;
                            break;
                        }
                    }
                }
            } else if (b2pos_of(kc, k) < seq.p) { // no dup after all
                *k |= DUP_BIT;
                --kc->reeval;
            }
            *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
            *k |= (uint64_t)!seq.t << ORIENT_SHFT | seq.p; // set new pos and strand
            get_next_nt_seq(s, seq);
            seq_next(seq);
        } while (++seq.p != p);
    }
}

static Hdr*
next_hdr(kct_t *kc, HK *hk, uint64_t *k)
{
    Hdr *h = kc->h + hk->hoffs;
    for (kc->bdit = h->bnd->end(); b2pos_of(*k) < (*--kc->bdit).s;)
        {}
    return h;
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
    uint64_t *kend = kc->kct + kc->kct_l - 1 - kc->last_uqct; // location after uniques, from last time
    HK *hk = kc->hk + kc->hk_l;
    Hdr *h = next_hdr(kc, --hk, kend);
    uint64_t *kctnext = kc->kct_next; // TODO: use this instead of pk below.
    std::set<uint64_t> pk;     // storage for unique key offsetss
    uint64_t *last = NULL;

    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    uint64_t *sk = kc->kct;
    for (uint64_t *k = kend; k != kc->kct - 1; --k) {

        if (IS_DUP(k))
            continue;

        pk.insert(*k);
        if (k - kc->kct < hk->koffs) {

            handle_range(kc, h, NULL, last, hk->koffs, &sk, pk);
            reached_boundary(kc, h, last, hk->koffs, pk);

            h = next_hdr(kc, --hk, k);
            last = NULL;
        }
        handle_range(kc, h, k, last, hk->koffs, &sk, pk);
        last = k;
    }
    handle_range(kc, h, NULL, last, hk->koffs, &sk, pk);
    reached_boundary(kc, h, last, hk->koffs, pk);

    place_uniques(kc, sk, kend, kc->s + h->s_s, pk);
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

    kc->kct_next = (uint64_t*)malloc(kc->kct_l * sizeof(uint64_t));
    if (kc->kct_next == NULL) {
        res = -ENOMEM;
        goto err;
    }
    for (kc->ext = res = 1; res > 0 && kc->ext != kc->readlength - KEY_WIDTH + 1;++kc->ext) {
        kc->iter = 0;
        do { // until no no more new uniques
            ext_uq_iter(kc);
        } while (kc->reeval > 0);
    }
    _buf_free(kc->kct_next);
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


