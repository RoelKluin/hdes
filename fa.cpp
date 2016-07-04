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
buf_grow_ks(kct_t *kc, Bnd &b, uint32_t **k1, uint32_t **k2)
{
    if ((kc->kct_l + 1) >= (1ul << kc->kct_m)) {
        unsigned ok1, ok3, ok4;

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
static inline uint32_t *
excise_one(kct_t *kc, Bnd &b, uint32_t *thisk, uint32_t *k)
{
    keyseq_t seq = {0};
    NB(k < kc->kct + kc->kct_l);
    buf_grow_ks(kc, b, &thisk, &k);
    NB(k != kc->kct + kc->kct_l);
    //*k &= ~DUP_BIT; // or keys that were moved after extension still have their dup bit set.
    kc->kct[kc->kct_l] = seq.p = *k;
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s, 0);
    *k ^= *k;//
    *contxt_idx = kc->kct_l++;
    ++b.moved;
    return thisk;//S;
}

static uint32_t *
excise(kct_t *kc, Bnd &b, uint32_t *thisk)
{
    for (uint32_t *k = b.tgtk + b.moved; k < thisk; ++k)
        thisk = excise_one(kc, b, thisk, k);

    return thisk;
}

//keys not in scope of unique
static inline void
move_uniq_one(kct_t *kc, Bnd &b, keyseq_t &seq, Mantra* bnd, uint32_t *contxt_idx)
{
    NB(*contxt_idx != NO_K);
    NB(*contxt_idx < kc->kct_l);
    uint32_t *k = kc->kct + *contxt_idx;
    if (k < b.tgtk) {
        //O; second+ occurance (may still be in scope)
        // after &&: do not set dupbit if previous key was within read scope (a cornercase)
        if ((~*k & DUP_BIT)/* && FIXME: issue in GE:8
                (k - kc->kct < b.fk || b2pos_of(*k) + b.ext < b2pos_of(seq.p))*/) {
            //P; no dup after all

            *k |= DUP_BIT;
            --kc->ct;
        }
    } else {

        //O; 1st occurance
        ++kc->ct;    // unique or decremented later.
        *k = seq.p; // set new pos and strand, unset dupbit
        if (b.tgtk != k) {
            NB(*k);
            *b.tgtk = *k;
            *k ^= *k;//
            *contxt_idx = b.tgtk - kc->kct;//K; moved up
        }
        ++b.tgtk;
        if (b.moved && b.tgtk[b.moved - 1])
            --b.moved;

        NB(b.tgtk <= kc->kct + kc->kct_l);
    }
}

static keyseq_t
move_uniq(kct_t *kc, Mantra* bnd, Bnd &b, C unsigned pend)
{
    C unsigned start = bnd->s;
    keyseq_t seq = {.p = start};
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, b.s);// already increments seq.p
    seq.p = b2pos_of(seq.p);
    NB(*contxt_idx < kc->kct_l);
    NB(*contxt_idx != NO_K);
    NB(seq.p <= pend);

    while (seq.p != pend) {
        move_uniq_one(kc, b, seq, bnd, contxt_idx);
        get_next_nt_seq(b.s, seq);
        contxt_idx = get_kct(kc, seq, 1);
        seq.p = b2pos_of(seq.p) + 2;
    }
    move_uniq_one(kc, b, seq, bnd, contxt_idx);

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
    t = kroundup32(t);
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
            EPR0(".");
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
                    k = excise(kc, b, k);
                    --k;
                } else {
                    //P; out of scope.
                    move_uniq(kc, bnd, b, end - 2);
                    if (end + 2 == bnd->e) { //prevent insertion & removal, + b.ext?
                        bnd->e = end;
                        ++k;
                        k = excise(kc, b, k);
                        break;
                    }
                    Mantra copy = *bnd;
                    copy.e = end;
                    buf_grow_add(kc->bnd, 1ul, 0, copy);
                    bnd->s = end + 2;
                    k = excise_one(kc, b, k, k);
                }
            }
            ++k;
        }

        if (bnd->s + b.ext >= bnd->e) {
            //P; in scope of start; excision (no bnd copy from b.obnd to kc->bnd)
            k = excise(kc, b, k);
        } else if (b2pos_of(*k) != bnd->e){
            //P; alt
            move_uniq(kc, bnd, b, bnd->e - 2);
            buf_grow_add(kc->bnd, 1ul, 0, *bnd);
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


