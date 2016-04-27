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
//#include <stack>
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"

static int
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
    return 0;
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

static inline unsigned
get_nextnt(kct_t C*C kc, uint64_t C p)
{
    return (kc->s[p>>2] >> ((p&3) << 1)) & 3;
}


// reset and remove non-uniqs, move uqs to start of kc->kct array
static inline int
extd_uq_by_p(kct_t *kc, pos_t p, pos_t pend, uint64_t ho, std::set<uint64_t> &pk)
{
    seq_t *ndxkct;
    int res;
    keyseq_t seq = {0};

    if (p + kc->ext + 1 >= pend && p + 1 != pend) {
        seq.p = p + ho - KEY_WIDTH;
        _build_key(kc, seq, p + ho);
        ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
        uint64_t*k = kc->kct + *ndxkct;
        reeval_keys_in_scope:

        for (;p != pend; ++p) {
            seq.t = get_nextnt(kc, p + ho);
            seq_next(seq);

            ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
            k = kc->kct + *ndxkct;

            // in next iteration: reeval whether this is still a dup.
            if (*k & DUP_BIT) {
                if (_B2POS_OF(kc, k) > p) {
                    ++kc->reeval;
                } else if (_B2POS_OF(kc, k) == p) {
                    // all remaining were excised this iteration
                    --kc->reeval;
                    // this position should be recognisable since it is the only with
                    // k before the adjacent uniq's k.
                    *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
                    *k |= ((uint64_t)!seq.t << ORIENT_SHFT) | p;
                    pk.insert(*k);
                }
            } else { // may be unique - next iter.
                ++kc->reeval;
            }
            if (IS_DBG_K(kc, k)) dbgndxkct = *ndxkct;
        }
        //ASSERT(IS_UQ(k), return _PRNT_SEQ_BY_POS(kc, p + ho - 1), Pfmt ", %lu", p, _B2POS_OF(kc, k));

        // need to determine kcts for inbetween. Same can occur twice, use end instead.
    }
    res = 0;
err:
    return res;
}

static int
reached_boundary(kct_t *kc, Hdr *h, uint64_t C*C last, std::set<uint64_t> &pk)
{
    C pos_t p = last ? B2POS_OF(*last) : (*kc->bdit).e;
    int res = 0;
    show_mantras(kc, h);
    if ((p > (*kc->bdit).s) && (p - kc->ext <= (*kc->bdit).s)) {
        if (last) {
            --kc->bdit;
        } else {
            kc->bdit = h->bnd->erase(kc->bdit);
        }
        _EVAL(extd_uq_by_p(kc, (*kc->bdit).s, p, h->s_s, pk));
    }
err:
    return res;    
}

static int swap_kct(kct_t *kc,  uint64_t *k1,  uint64_t *k2, seq_t *ndxkct2, uint64_t ho)
{
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {0};
    pos_t p = _B2POS_OF(kc, k1);
    seq.p = p + ho - KEY_WIDTH;
    _build_key(kc, seq, p + ho);
    seq_t *ndxkct1 = _get_kct(kc, seq, seq.t, return -EFAULT);

    // ndxkct points to reordered kcts
    uint64_t t = *ndxkct1;
    *ndxkct1 = *ndxkct2;
    *ndxkct2 = t;
    
    // also swap info at kcts
    t = *k1;
    *k1 = *k2;
    *k2 = t;
    return 0;
}

// XXX XXX hk need insertions.
static int
place_uniques(kct_t *kc, uint64_t *sk, uint64_t* &kend, uint64_t ho, std::set<uint64_t> &pk)
{
    pos_t p;
    int res;
    HK *hk = kc->hk;
    ASSERT(0, return -EFAULT, "FIXME: hk need insertions. koffs need to be updated");

    for (std::set<uint64_t>::iterator it = pk.begin(); it != pk.end(); ++it) {
        uint64_t kv = *it;
        if (~kv & DUP_BIT) {
            ++kc->uqct;
            p = B2POS_OF(kv);
            keyseq_t seq = {0};
            seq.p = p + ho - KEY_WIDTH;
            _build_key(kc, seq, p + ho);
            seq_t *ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
            // If smaller than sk, an ambiguous key was apparently unique. These are handled in
            // next iteration
            ASSERT(kc->kct + *ndxkct >= sk, return -EFAULT);

            _EVAL(swap_kct(kc, kc->kct + *ndxkct, kend--, ndxkct, ho));
        }
    }

err:
    return res;
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
    if (prevk) {
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


static int
handle_range(kct_t *kc, Hdr *h, uint64_t C*C prevk, uint64_t C*C thisk,
        C uint32_t koffs, uint64_t** sk, std::set<uint64_t> &pk)
{
    int res;
    seq_t *ndxkct;
    keyseq_t seq = {0};
    uint64_t ho = h->s_s;
    C pos_t b2start = (*kc->bdit).s;
    C pos_t b2end = (*kc->bdit).e;
    pos_t prev = prevk ? B2POS_OF(*prevk) : b2start + KEY_WIDTH - 1;
    C pos_t p = thisk ? B2POS_OF(*thisk) : b2end;

    if (p - prev <= kc->ext + 1) { // a 2nd uq
        uq_bounded(kc, h, prevk, thisk, prev, p);

        // what if all remaining occurances happen between two boundaries?
        if (prev - p > 1)
            _EVAL(extd_uq_by_p(kc, prev, p, ho, pk));
        prev = p;
        seq.p = prev + ho - KEY_WIDTH;
        _build_key(kc, seq, prev + ho);
        ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);

    } else { // first occurance.

        // 1: from prev till now, add position and dup reevaluation
        seq.p = ho + prev + 1 - KEY_WIDTH;
        _build_key(kc, seq, ho + prev + 1);
        do {
            ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
            uint64_t *k = kc->kct + *ndxkct;

            if (*k & DUP_BIT) {
                if (k - kc->kct < koffs || _B2POS_OF(kc, k) + ho < seq.p) {
uq_1st_occurance:   *k &= ~DUP_BIT; // unset for 1st occurance
//                    _EVAL(swap_kct(kc, *sk, k, ndxkct, ho)); //XXX XXX: valid ???
                    ++kc->reeval;
//                    k = (*sk)++;
                } else if (B2POS_OF(*k) == seq.p - ho) {
                    // all but last positions were excised, unique after all!
                    uq_bounded(kc, h, prevk, thisk, prev, seq.p - ho);
                    pk.insert(*k);
                    unsigned at = seq.p - ho;
                    if (at < prev + kc->ext + 1) {
                        if (at > prev + 1)
                            _EVAL(extd_uq_by_p(kc, prev, at, ho, pk));
                        prev = at;
                        if (p - prev < kc->ext + 1) {
                            //rollback entirely
                            uq_bounded(kc, h, prevk, NULL, prev, p);
                            seq.p = p + ho - KEY_WIDTH;
                            _build_key(kc, seq, p + ho);
                            ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
                            if (p - prev > 1)
                                _EVAL(extd_uq_by_p(kc, prev, p, ho, pk));
                            prev = p;
                            break;
                        }
                    }
                }
            } else if (_B2POS_OF(kc, k) < seq.p - ho) { // no dup after all
                *k |= DUP_BIT;
                --kc->reeval;
            }
            *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
            *k |= (uint64_t)!seq.t << ORIENT_SHFT | (seq.p - ho); // set new pos and strand
            seq.t = get_nextnt(kc, seq.p);
            seq_next(seq);
        } while (++seq.p != p + ho);
    }
    res = 0;
err:
    if (res < -1) {
        EPR("prev:" Pfmt ", p:" Pfmt ", b2end:" Pfmt, prev, p, b2end);
    }
    return res;
}

static Hdr*
next_hdr(kct_t *kc, HK *hk, uint64_t *k)
{
    Hdr *h = kc->h + hk->hoffs;
    for (kc->bdit = h->bnd->end(); B2POS_OF(*k) < (*--kc->bdit).s;)
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
static int
ext_uq_iter(kct_t *kc)
{
    int res;

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

            _EVAL(handle_range(kc, h, NULL, last, hk->koffs, &sk, pk));
            _EVAL(reached_boundary(kc, h, last, pk));

            h = next_hdr(kc, --hk, k);
            last = NULL;
        }
        _EVAL(handle_range(kc, h, k, last, hk->koffs, &sk, pk));
        last = k;
    }
    _EVAL(handle_range(kc, h, NULL, last, hk->koffs, &sk, pk));
    _EVAL(reached_boundary(kc, h, last, pk));

    _EVAL(place_uniques(kc, sk, kend, h->s_s, pk));
    //ASSERT(0, return -EFAULT, "(NOT!;) end of loop reached! (TODO: swapping..)");
    // FIXME: iteration over swapped instead of uniques
    //

    // TODO: some uniques no longer occur: all are excised.
    kc->last_uqct = kc->uqct;

    EPQ(dbg > 0, "observed %u uniques in iteration %u, extension %u, %u to be reevaluated\n",
            kc->uqct, ++kc->iter, kc->ext, kc->reeval);

    res = kc->reeval;
err:
    return res;
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
            res = ext_uq_iter(kc);
        } while (res > 0);
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
        ASSERT(fh[2].name != NULL, return -EFAULT);
        fh[0].name = file;

        len = strlen(fh[2].name);
        strncpy(fh[0].name, fh[2].name, ++len);
        _ACTION0(reopen(fh, ext[0], ext[5]), "");
    }
    ASSERT(len < 256, return -EFAULT, "filename too long: %s", fh[0].name);
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


