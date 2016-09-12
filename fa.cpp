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
free_kc(Key_t *kc)
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
print_kct(Key_t *kc, Mantra* at, Ext_t* e, uint32_t* tk)
{
    unsigned i = 0;
    EPR("");
    Hdr* h = kc->h;
    uint8_t *s = kc->s;
    Mantra* bnd = kc->bnd_l ? kc->bnd : at;
    uint32_t* hkoffs = kc->hkoffs;

    for (uint32_t*k = kc->kct; k - kc->kct < kc->kct_l;) {
        if (bnd != e->obnd + e->obnd_l) {
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
            char c = k!=tk?(k!=e->tgtk?' ':'t'):'k';
            EPR0("%c>%s:%u\t%u\t%c", c, kc->id + h->ido, i, *hkoffs, *k & DUP_BIT?'*':' ');
            if (*k != 0)
                print_posseq(s, *k);
            else
                EPR("(removed)");

            ++k;
        }
        if (bnd != e->obnd + e->obnd_l && ++bnd == kc->bnd + kc->bnd_l)
            bnd = at;
    }
}

// if kct grows, pointers *thisk, k, e->tgtk become invalid
static void
buf_grow_ks(Key_t *kc, Ext_t* e, unsigned add, uint32_t **k)
{
    while ((kc->kct_l + add) >= (1ul << kc->kct_m)) {
        EPR("buf_grow_ks grew");
        unsigned ok1, ok2;

        ok1 = *k - kc->kct;
        ok2 = e->tgtk - kc->kct;
        uint32_t *t = (uint32_t *)realloc(kc->kct, sizeof(uint32_t) << ++kc->kct_m);
        if_ever (t == NULL) {
            EPR("buf_grow_ks failed");
            raise(SIGTRAP);
        }
        kc->kct = t;
        *k = t + ok1;
        e->tgtk = t + ok2;
    }
}

/*
 * e->obnd contains ranges per contig for sequence not yet be mappable, but may become so.
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
k_compression(Key_t *kc, Ext_t* e, uint32_t *hkoffs, uint32_t *k)
{

    NB(hkoffs <= kc->hkoffs + kc->h_l);
    while (hkoffs != kc->hkoffs + kc->h_l)
        *hkoffs++ = e->tgtk - kc->kct;

    Hdr* h = kc->h;
    e->s = kc->s;

    while (k - kc->kct != kc->kct_l) {
        NB(hkoffs < kc->hkoffs + kc->hkoffs_l);
        while (k - kc->kct == *hkoffs) {
            if (h - kc->h != kc->h_l - 1) {
                e->s += h++->len;
            } else {
                h = kc->h;
                e->s = kc->s;
            }
            *hkoffs++ = e->tgtk - kc->kct;
        }
        if (*k) {
            if (k != e->tgtk) {
                keyseq_t seq = {.p = b2pos_of(*k) };
                kc->contxt_idx[build_ndx_kct(kc, seq, e->s, 0)] = e->tgtk - kc->kct;
                *e->tgtk = *k;
                *k ^= *k;//
            }
            ++e->tgtk;
        }
        ++k;
    }
    kc->kct_l = e->tgtk - kc->kct;
    while (hkoffs != kc->hkoffs + kc->hkoffs_l)
        *hkoffs++ = kc->kct_l;
}

// keys between two uniqs.
static inline void
excise_one(Key_t *kc, Ext_t* e, uint32_t *q)
{
    keyseq_t seq = {0};
    NB(q < kc->kct + kc->kct_l);
    NB(q != kc->kct + kc->kct_l);
    //*q &= ~DUP_BIT; // or keys that were moved after extension still have their dup bit set.
    kc->kct[kc->kct_l] = seq.p = *q;
    uint32_t *contxt_idx = kc->contxt_idx + build_ndx_kct(kc, seq, e->s, 0);
    *q ^= *q;//
    *contxt_idx = kc->kct_l++;
    ++e->moved;//S;
}

// keys not in scope of unique. hot.
// * can we use a Ext_t.moved and move all keys if none are uniq?
static void
move_uniq(Key_t *kc, Ext_t* e, C unsigned start, C unsigned pend)
{
    unsigned dna = 0, rc = 0, t, key_complete = start - 2, p = start - NT_WIDTH;
    while (p != key_complete) { // build key
        t = (e->s[p>>3] >> (p&6)) & 3;
        rc = ((rc << 2) & KEYNT_MASK) | (t ^ 2);
        dna = t << KEYNT_TOP | dna >> 2;
        p += 2;
    }
    do {
        t = (e->s[p>>3] >> (p&6)) & 3;
        rc = ((rc << 2) & KEYNT_MASK) | (t ^ 2);
        dna = t << KEYNT_TOP | dna >> 2;

        // build complement independent index (orientation bit extracted)
        unsigned dev = dna ^ rc;
        t = dev & -dev;            /* isolate deviant bit */
        t |= !t;                   /* for palindromes: have to set one. */
        unsigned ori = !(t & dna); /* orientation: was devbit set? */
        t = dna ^ (-ori & dev);    /* index is based on dna or rc dependent on orientation */
        unsigned top = t & KEYNT_BUFSZ;  /* if the top bit is set.. */
        t ^= (-!!top) & SNDX_TRUNC_MASK; /* ..flip all bits to shorten index by one */

        p += 2;
        uint32_t *k = kc->kct + kc->contxt_idx[t];
        if (k > e->tgtk) {
            //O; 1st occurance
            ++kc->ct;           // unique or decremented later.
            *k = p | ori;       // set new pos and strand, unset dupbit

            *e->tgtk = *k;
            *k ^= *k;//
            kc->contxt_idx[t] = e->tgtk - kc->kct;
            e->moved -= (e->tgtk[e->moved] != 0);
            ++e->tgtk;
        } else if (k < e->tgtk) {
            //O; second+ occurance (may still be in scope)
            // after &&: do not set dupbit if previous key was within read scope (a cornercase)
            t = (~*k & DUP_BIT) && (k - kc->kct < e->fk || b2pos_of(*k) + kc->ext < p);
            //P; no dup after all if t is set.

            *k |= -t & DUP_BIT;
            kc->ct -= t;
        } else {

            //O; also 1st occurance
            ++kc->ct;     // unique or decremented later.
            *k = p | ori; // set new pos and strand, unset dupbit

            ++e->tgtk;
        }
    } while (p != pend);
}

/*
 * Called for extensions 1+, primary uniques were already determined upon reading the fasta.
 * secondary uniques arise when all but one key are already covered by primary unique keys.
 *
 * 1) iterate over as of yet ambiguous keys, determine which ones have become unique.
 * 2) shrink remaining ambiguous regions (mantra) in scope of unique keys.
 *
 * keys are kept ordered. Ambiguous keys are kept ordered upon first occurance. Unique keys are
 * isolated and ordered on extension, contig and thirdly on position.
 *
 */
static void
ext_uq_iter(Key_t *kc, Ext_t* e)
{
    uint32_t *k = kc->kct;
    uint32_t *hkoffs = kc->hkoffs;
    Mantra *bnd = kc->bnd;
    Hdr *h = kc->h;
    unsigned skctl = kc->kct_l;

    kc->bnd = e->obnd; // buffer will be reused
    e->obnd = bnd; // swap e->obnd and kc->bnd.
    e->obnd_l = kc->bnd_l;
    e->obnd_m = kc->bnd_m;

    unsigned t = e->obnd_l;
    kroundup32(t);
    kc->bnd_m = __builtin_ctz(t) + 1;
    kc->bnd_l = 0;

    buf_realloc(kc->bnd, kc->bnd_m);
    //bnd = e->obnd; (already)

    e->tgtk = k;
    e->s = kc->s;
    e->fk = k - kc->kct;
    e->moved = 0;//B; initial state

    // TODO: skip k up to first promising in previous iteration.
    do {
        NB(h - kc->h <= bnd->ho);
        while (h - kc->h != bnd->ho) {
            //~ also update header
            *hkoffs++ = e->tgtk - kc->kct;
            t = hkoffs - kc->hkoffs;
            buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
            hkoffs = kc->hkoffs + t;
            e->s += h++->len;
            EPR0(".");
            e->fk = k - kc->kct;
        }
        NB(hkoffs == kc->hkoffs + bnd->ho);

        //k - kc->kct <= *hkoffs: may be untrue after key excision.
        //B; next region

        while ((k - kc->kct) < *hkoffs && b2pos_of(*k) < bnd->e) {

            if (IS_UQ(k)) { //~ uniq
                unsigned end = b2pos_of(*k);
                if (bnd->s + kc->ext < end) {
                    //P; out of scope.
                    move_uniq(kc, e, bnd->s, end - 2);
                    Mantra copy = *bnd;
                    copy.e = end;
                    buf_grow_add(kc->bnd, 1ul, 0, copy);
                    bnd->s = end + 2;
                    buf_grow_ks(kc, e, 1, &k);
                    excise_one(kc, e, k);
                } else {
                    //P; in scope of start; excision
                    bnd->s = end + 2;
                    ++k;
                    buf_grow_ks(kc, e, (k - e->tgtk) - e->moved, &k);
                    for (uint32_t *q = e->tgtk + e->moved; q < k; ++q)
                        excise_one(kc, e, q);
                    --k;
                }
            }
            ++k;
        }

        if (bnd->s + kc->ext >= bnd->e) {
            //P; in scope of start; excision (no bnd copy from e->obnd to kc->bnd)
            buf_grow_ks(kc, e, (k - e->tgtk) - e->moved, &k);
            for (uint32_t *q = e->tgtk + e->moved; q < k; ++q)
                excise_one(kc, e, q);
        } else {
            //P; alt
            move_uniq(kc, e, bnd->s, bnd->e - 2);
            buf_grow_add(kc->bnd, 1ul, 0, *bnd);//K;
        }
    } while (++bnd != e->obnd + e->obnd_l);

    t = hkoffs - kc->hkoffs;
    while (h - kc->h != kc->h_l) {
        buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
        ++h;
    }
    k_compression(kc, e, kc->hkoffs + t, k);//K;

    NB(skctl == kc->kct_l, "skctl(%u) != kc->kct_l(%u)", skctl, kc->kct_l);//B; final state
}

static int
extd_uniqbnd(Key_t *kc)
{
    int res;
    Ext_t e = {0};
    // FIXME: don't stop at a certain readlength, stop if there are no more dups or
    // no more that can become one
    unsigned end = (kc->readlength - KEY_WIDTH + 1) << 1;
    kc->ext_iter = buf_init(kc->ext_iter, 1);
    for (kc->ext = 0; kc->ext != end; kc->ext += 2) {
        unsigned iter = 0;
        if (kc->hkoffs[kc->h_l-1] != 0) { // or all keys were already finished.
            do { // until no no more new uniques
                kc->ct = 0;
                ext_uq_iter(kc, &e);
                EPR("observed %u potential in iteration %u, extension %u\n",
                    kc->ct, ++iter, kc->ext >> 1);
            } while (kc->ct > 0);
        }
        EPR("----[ end of extension %u ]-------", kc->ext >> 1);
        buf_grow_add(kc->ext_iter, 1ul, 0, iter);
    }
    res = 0;
err:
    free(e.obnd);
    return res;
}

int
fa_index(struct gzfh_t *fh, uint64_t optm, unsigned readlength)
{
    int len, res = -ENOMEM;
    char file[512];
    Key_t kc = {0};
    kc.readlength = readlength;
    unsigned mode;

    C char *ext[4] = {".fa", ".2b", ".kc", ".uq"}; // .kc is a temp file

    if (fh[0].name) {
        len = strlen(fh[0].name);
    } else {
        NB(fh[2].name != NULL);
        fh[0].name = file;

        len = strlen(fh[2].name);
        strncpy(fh[0].name, fh[2].name, ++len);
        _ACTION0(reopen(fh, ext[0], ext[3]), ""); // .fa => .uq
    }
    NB(len < 256, "filename too long: %s", fh[0].name);
    if (fh[1].name == NULL)
        fh[1].name = &file[len];

    strncpy(fh[1].name, fh[0].name, len);

    kc.contxt_idx = buf_init_arr(kc.contxt_idx, KEYNT_BUFSZ_SHFT);
    // first check whether unique keys are there already.
    if (fh[0].fp) {
        mode = 2;
    } else {
        _ACTION0(reopen(fh, ext[3], ext[2]), "No %s trying to open %s instead", ext[3], ext[2]); // uq => kc
        mode = fh[0].fp != NULL;
    }
    _ACTION0(reopen(fh + 1, ext[3], ext[1]), "%s does%s exist", ext[1], fh[1].fp ? "" : " not");

    if (mode && fh[0].fp && fh[1].fp) {
         _ACTION(load_seqb2(fh + 1, &kc), "loading twobit sequence file");
         _ACTION(load_kc(fh, &kc), "loading keycounts file");
    } else {
        bool found = false;
        for (int i=0; i != 2; ++i) {
            if (fh[i].fp) {
                if (~optm & amopt('f')) {
                    found = true;
                    EPR("%s file present, refusing to overwrite.", fh[i].name);
                }
                rclose(fh +i);
            }
        }
        if (optm & amopt('f')) {
            mode = 0;
        } else if (mode) {
            EPR("%s file present, refusing to overwrite.", ext[3+mode]);
            goto err;
        }
        if (found) goto err;
        EPR("starting from scratch.");
        _ACTION(fa_read(fh, &kc), "reading fasta");
    }
    if (mode < 2) {
        _ACTION(reopen(fh, ext[2], ext[3]), "");
        _ACTION(extd_uniqbnd(&kc), "continuing kc");
        _ACTION(save_kc(fh, &kc), "writing unique keycounts file");
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


