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
    buf_free(kc->hk);
    buf_free(kc->ext);
    buf_free(kc->contxt_idx);
}

static void
print_kct(kct_t *kc, Bnd &b, uint32_t* tk)
{
    unsigned i = 0;
    unsigned offs = 0;
    uint8_t *s = kc->s;
    HK *hk;
    EPR("");
    for (hk = kc->hk; hk != kc->hk + kc->hk_l; ++hk) {
        while (i < hk->koffs) {
            uint32_t* k = kc->kct + i;
            char c = k!=tk?(k!=b.prev?(k!=b.sk?' ':'s'):'p'):'k';
            EPR0("%u%c:\t0x%x\t", i, c, kc->kct[i]);
            if (kc->kct[i] != 0)
                print_posseq(s, b2pos_of(kc->kct[i]));
            else
                EPR("(removed)");
            ++i;
        }
        s += (hk->len >> 2) + !!(hk->len & 3);
    }
    uint64_t* ext = kc->ext;
    while (i < kc->kct_l) {
        s = kc->s;
        for (hk = kc->hk; hk != kc->hk + kc->hk_l; ++hk) {
            unsigned j = ext - kc->ext < kc->ext_l ? *ext++ : ~0u;
            while (j) {
                uint32_t* k = kc->kct + i;
                char c = k!=tk?(k!=b.prev?(k!=b.sk?' ':'s'):'p'):'k';
                EPR0("~%u%c:\t0x%x\t", i, c, kc->kct[i]);
                if (kc->kct[i] != 0) {
                    --j;
                    print_posseq(s, b2pos_of(kc->kct[i]));
                } else {
                    EPR("(removed)");
                }
                if (++i == kc->kct_l)
                    return;
            }
            s += (hk->len >> 2) + !!(hk->len & 3);
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
shrink_mantra(kct_t *kc, Bnd &b, uint32_t C*C thisk, C uint32_t prev, C uint32_t p)
{
    NB(b.it != kc->bnd->end());
    if (b.prev && prev != (*b.it).s) {     // not at mantra start
        C uint32_t end = (*b.it).e;
        C uint32_t ke = (*b.it).ke;
        (*b.it).e = prev;
        (*b.it).ke = b.sk - kc->kct;
        if (thisk - kc->kct != (*b.it).ke) {         // not at mantra end either
            kc->bnd->insert(b.it, *b.it);  // copy of current
            (*b.it).s = p;                 // mantra became two smaller ranges.
            (*b.it).e = end;
            (*b.it).ke = ke;
        } // or end got shifted
    } else {
        // shift start or if at end, force erase (entirely mapable).
        (*b.it).s = (thisk - kc->kct != (*b.it).ke) ? p : (*b.it).e; // or entire region became mapable, force erase.
    }
    NB((*b.it).s <= (*b.it).e);
    if ((*b.it).s == (*b.it).e) {
        b.it = kc->bnd->erase(b.it);   // nothing left
        b.prev = NULL;
    }//GDB:mantra1
}

static void
process_mantra(kct_t *kc, Bnd &b, uint32_t *C thisk)
{
    C uint32_t prev = _prev_or_bnd_start(b);
    C uint32_t pend = (thisk - kc->kct != (*b.it).ke) ? b2pos_of(*thisk) - 1 : (*b.it).e;

    uint32_t *contxt_idx;
    keyseq_t seq = {0};

    if (in_scope(kc, prev, pend)) { // a 2nd uniq
        shrink_mantra(kc, b, thisk, prev, pend);
        // The distance between b.sk and k may have grown by excised 1st kcts (beside uniq).
        if (thisk - b.sk > b.moved) {
            EPR("%u were moved during excision", (thisk - b.sk) - b.moved + 1);

            for (uint32_t *k = b.sk + b.moved; k <= thisk; ++k) {
                NB(k < kc->kct + kc->kct_l);
                seq.p = b2pos_of(kc, k);
                contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s);
                buf_grow(kc->kct, 1, 0);
                kc->kct[kc->kct_l] = *k;
                *k ^= *k;//
                *contxt_idx = kc->kct_l++;//GDB:1
            }
            b.moved = b.sk - thisk + 1;
            b.prev = kc->kct + *contxt_idx;
        } else {
            b.prev = thisk;
        }
        return;
    }
    // prev to pend are uniques, not in scope. between uniqs are
    // from prev to p, add position if pending and reevaluate dupbit
    seq.p = prev + 1;
    contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s); // already increments seq.p
    NB(*contxt_idx != NO_KCT);
EPR("kept %u-%u", prev, pend);

    for (;;) {
        uint32_t *k = kc->kct + *contxt_idx;
print_seq(&seq);

        if (k < b.sk) { // XXX requires multi-contig uniqs and excised in b.sk .. k
            // second occurance

            if (~*k & DUP_BIT) {
                EPR("// no dup after all");

                *k |= DUP_BIT;
                --kc->uqct;
            } else {
                EPR("dup bit was already set");
            }
        } else {
            EPR("// 1st occurance");
            ++kc->uqct;    // unique or decremented later.

            if (k - kc->kct >= (*b.it).ke)
                EPR("// a position is pending for %u'th (!= %u), first was excised", seq.p, b2pos_of(kc, k));

            *k = seq.p << 1 | (seq.t != 0); // set new pos and strand, unset dupbit
            if (b.sk != k) {
                *b.sk = *k;
                *k ^= *k;//
                *contxt_idx = b.sk - kc->kct;//GDB:move
                EPR("moved");
            }
            ++b.sk;
            NB(b.sk <= kc->kct + kc->kct_l);
        }
        if (seq.p == pend)
            break;

        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq);
        ++seq.p;
    }
    if (thisk - kc->kct != (*b.it).ke) {
        //NB((thisk - b.sk) - b.moved == 0, "%u - %u", (b.sk - thisk), b.moved);
        EPR("only one uniq isolated from mantra");
        ++b.moved; //namely last uniq.
        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq);
        buf_grow(kc->kct, 1, 0);
        uint32_t *k = kc->kct + kc->kct_l;
        if (pend != (*b.it).e) {
            kc->bnd->insert(b.it, *b.it);  // copy of current
            (*--b.it).e = pend;
            (*b.it).ke = b.sk - kc->kct;
            (*++b.it).s = pend + 1;
            //mantra4
        }
        NB(k != thisk);
        *k = *thisk;
        *thisk ^= *thisk;//
        *contxt_idx = kc->kct_l++;//GDB:2
        b.prev = k;
    } else {
        b.prev = thisk;
    }
}

static void
reached_boundary(kct_t *kc, Bnd &b)
{
    C uint32_t prev = _prev_or_bnd_start(b);

    // if at start and entirely mapable: remove mantra.
    NB(b.it != kc->bnd->end());
    if (b.prev == NULL && in_scope(kc, prev, (*b.it).e)) {
EPR("mantra removed");//GDB:mantra2
        b.it = kc->bnd->erase(b.it);
    } else {
EPR("next bnd");//GDB:mantra3
        (*b.it).ke = b.sk - kc->kct;
        ++b.it;
    }
}

static inline void
skip_mantra(kct_t *kc, Bnd &b, uint32_t *k)
{
    process_mantra(kc, b, k);
    if (b.it != kc->bnd->end())
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
    HK *hk;
    uint32_t *k = kc->kct; // location after uniques, from last time
    Bnd b = {
        .sk = k,
        .s = kc->s,
        .prev = NULL,
        .moved = 0,
        .it = kc->bnd->begin()
    };

    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    for (hk = kc->hk; hk != kc->hk + kc->hk_l; ++hk) {

        if (k < kc->kct + hk->koffs) {
EPR("%u k's left on bnd", hk->koffs - (k - kc->kct));

            while (k < kc->kct + hk->koffs) {

                while (k - kc->kct < (*b.it).ke) {
EPR("%u", k - kc->kct);

                    if (IS_UQ(k)) { // they are handled during excision (or postponed) - process_mantra()
EPR("uniq");
print_posseq(b.s, b2pos_of(*k));
                        process_mantra(kc, b, k); //
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

        buf_grow_add(kc->ext, 1ul, 0, uq_and_1stexcised);
        EPQ(uq_and_1stexcised, "total %u uniq", uq_and_1stexcised);
        hk->koffs = b.sk - kc->kct;
        //kc->uqct = uq_and_1stexcised;
    }
out:
    kc->uqct = k - b.sk;
    kc->last_uqct = kc->uqct;
    b.s = kc->s;
    uint64_t* ext = kc->ext;
    // put uniqs and 1st excised back in array (recompression)
    for (hk = kc->hk; hk != kc->hk + kc->hk_l; ++hk) {
        int j = *ext++;
        while (j--) {
            if (*k != 0) {
                keyseq_t seq = {.p = b2pos_of(kc, k) };
                uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s);
                *contxt_idx = b.sk - kc->kct;
                *b.sk++ = *k;
                *k ^= *k;//
            }
            ++k;
        }
        b.s += (hk->len >> 2) + !!(hk->len & 3);
    }
    kc->kct_l = b.sk - kc->kct;
    return 0;
}

static int
extd_uniqbnd(kct_t *kc, struct gzfh_t *fhout)
{
    int res = -ENOMEM;
    kc->ext = buf_init(kc->ext, 1);
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


