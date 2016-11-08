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

extern void
print_kct(Key_t *kc, Mantra* at, Iter_t* e, uint32_t* tk)
{
    unsigned extension = 0;
    EPR("");
    Hdr* h = kc->h;
    uint8_t *s = kc->s;
    Mantra* bnd = kc->bnd_l ? kc->bnd : at;
    uint32_t* hkoffs = kc->hkoffs;

    for (uint32_t*k = kc->kct; K_OFFS(kc, k) < kc->kct_l;) {
        // print correctly, using old boundaries, if print_kct is called during (ext_uq) iter.
        if (bnd != e->obnd + e->obnd_l) {
            while (h - kc->h != bnd->ho)
                s += h++->len;
            hkoffs = kc->hkoffs + bnd->ho;
        } else {
            if (h - kc->h != kc->h_l - 1) {
                s += h++->len;
            } else {
                ++extension;
                h = kc->h;
                s = kc->s;
            }
            if (++hkoffs == kc->hkoffs + kc->hkoffs_l)
                hkoffs = &kc->kct_l;
        }

        while (K_OFFS(kc, k) < *hkoffs) {
            char c = k!=tk?(k!=e->tgtk?' ':'t'):'k';
            EPR0("%c>%s:%u\t%u\t%c", c, kc->id + h->ido, extension,
                    *hkoffs, *k & DUP_BIT?'*':' ');
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
buf_grow_ks(Key_t *kc, Iter_t* e, unsigned add, uint32_t **k)
{
    while ((kc->kct_l + add) >= (1ul << kc->kct_m)) {
        EPR("buf_grow_ks grew");
        unsigned ok1, ok2;

        ok1 = K_OFFS(kc, *k);
        ok2 = K_OFFS(kc, e->tgtk);
        kc->kct = (uint32_t *)realloc(kc->kct, sizeof(uint32_t) << ++kc->kct_m);
        NB (kc->kct != NULL);
        *k = kc->kct + ok1;
        e->tgtk = kc->kct + ok2;
    }
}

/*
 * e->obnd contains ranges per contig for sequence not yet be mappable, but may become so.
 * Initially there are one or more dependent on the presence and number of N-stretches;
 * initialized in key_init.cpp, But as uniques within read scope of one another cover
 * regions, this remaining `mantra' shrinks.
 *
 * mantras start at the chr start, the position after a N-stretch or after a region 'scoped'
 * by uniques, the first non-unique after. The mantra ends at chr end, stretch or first non
 * unique before a region scoped by uniques.
 *
 * while the no. regions can grow, the total no. sequence decreases.
 */

/*
 * During an iteration over the genome unique keys were excised and added to the end. The
 * causes a sparse k-mer array. hot!
 */
static void
denser_k(Key_t *kc, Iter_t* e, uint32_t *hkoffs, uint32_t *k)
{
    Hdr* h = kc->h;
    e->s = kc->s;


    // The k's of which 1st pos lied between uniques, but then appeared still recurrent,
    // on genome may have caused holes in the k array (kc->kct) these holes are undone here.
    while (K_OFFS(kc, k) != kc->kct_l) {
        NB(hkoffs < kc->hkoffs + kc->hkoffs_l);
        while (K_OFFS(kc, k) == *hkoffs) { // update k offsets
            if (h - kc->h != kc->h_l - 1) {
                e->s += h++->len;
            } else {
                h = kc->h;
                e->s = kc->s;
            }
            *hkoffs++ = K_OFFS(kc, e->tgtk);
        }
        if (*k) { // not a 'hole'
            if (k != e->tgtk) {
                keyseq_t seq = {.p = b2pos_of(*k) };
                kc->contxt_idx[build_ndx_kct(kc, seq, e->s, 0)] = K_OFFS(kc, e->tgtk);
                *e->tgtk = *k;
                *k ^= *k;//
            }
            ++e->tgtk;
        }
        ++k;
    }
    kc->kct_l = K_OFFS(kc, e->tgtk);
    while (hkoffs != kc->hkoffs + kc->hkoffs_l)
        *hkoffs++ = kc->kct_l;
}

// keys between two uniqs.
static inline void
excise_one(Key_t *kc, Iter_t* e, uint32_t *q)
{
    keyseq_t seq = {0};
    NB(K_OFFS(kc, q) < kc->kct_l);
    //*q &= ~DUP_BIT; // or keys that were moved after extension still have their dup bit set.
    kc->kct[kc->kct_l] = seq.p = *q;
    kc->contxt_idx[build_ndx_kct(kc, seq, e->s, 0)] = kc->kct_l++;
    *q ^= *q;//
    ++e->moved;//S;
}

// keys not in scope of unique. hot.
// * can we use a Iter_t.moved and move all keys if none are uniq?
static void
move_uniq(Key_t *kc, Iter_t* e, C unsigned start, C unsigned pend)
{
    unsigned dna = 0, rc = 0;

    for (unsigned p = start - NT_WIDTH; p != pend;) {

        unsigned t = (e->s[p>>3] >> (p&6)) & 3;
        rc = ((rc << 2) & KEYNT_MASK) | (t ^ 2);
        dna = t << KEYNT_TOP | dna >> 2;
        p += 2;

        // once key is complete
        if (p >= start) {

            // build complement independent index (orientation bit extracted)
            rc ^= dna;                  /* becomes deviant */
            t = rc & -rc;               /* isolate deviant bit */
            unsigned ori = !((t | !t) & dna);  /* orientation: 1st devbit or 1st bit for palidromes */
            t = dna ^ (-ori & rc);      /* index is based on dna or rc dependent on orientation */
            rc ^= dna;
            unsigned top = t & KEYNT_BUFSZ;  /* if the top bit is set.. */
            t ^= (-!!top) & SNDX_TRUNC_MASK; /* ..flip all bits to shorten index by one */

            uint32_t *k = kc->kct + kc->contxt_idx[t];
            if (k >= e->tgtk) {
                //O; 1st occurance
                ++kc->ct;           // unique or decremented later.
                *k = p | ori;       // set new pos and strand, unset dupbit

                if (k > e->tgtk) { // one earlier, this iteration, was in scope and excised.
                    *e->tgtk = *k;
                    *k ^= *k;//
                    kc->contxt_idx[t] = K_OFFS(kc, e->tgtk);
                    if (e->tgtk[e->moved])
                        --e->moved;
                }
                ++e->tgtk;
            } else {
                // O; second+ occurance, but can still be both in same scope (:after &&)!
                t = (~*k & DUP_BIT) && (K_OFFS(kc, k) < e->fk || b2pos_of(*k) + kc->ext < p);
                //P; no dup after all if t is set. TODO: move in this case!

                *k |= -t & DUP_BIT;
                kc->ct -= t;
            }
        }
    }
}

/*
 * Called for extensions 1+, primary uniques were already determined when the fasta was read.
 * Secondary uniques arise when all but one key are covered by primary unique keys in scope
 * of one another.
 *
 * 1) Iterate over as of yet ambiguous keys, determine which ones have become unique, this iter.
 * 2) shrink the region hence covered (decrease mantra).
 *
 * keys are kept ordered. Ambiguous keys are kept ordered upon first occurance. Unique keys are
 * isolated and ordered on extension, contig and thirdly on position.
 */
static void
ext_uq_iter(Key_t *kc, Iter_t* e)
{
    uint32_t *k = kc->kct;

    uint32_t *hkoffs = kc->hkoffs;
    Hdr *h = kc->h;

    Mantra *bnd = e->obnd;

    e->tgtk = k;
    e->s = kc->s;
    e->fk = K_OFFS(kc, k);
    e->moved = 0;//B; initial state

    kc->ct = 0;

    // TODO: skip k up to first promising in previous iteration.
    do {

        NB(h - kc->h <= bnd->ho);
        while (h - kc->h != bnd->ho) {
            //~ also update header
            hkoffs = update_append_hkoffs(kc, e, hkoffs);
            e->s += h++->len;
            e->fk = K_OFFS(kc, k);
            EPR0(".");
        }
        NB(hkoffs == kc->hkoffs + bnd->ho);

        //K_OFFS(kc, k) <= *hkoffs: may be untrue after key excision.
        //B; next region

        while (K_OFFS(kc, k) < *hkoffs && b2pos_of(*k) < bnd->e) {

            if (IS_UQ(k)) { //~ uniq
                unsigned end = b2pos_of(*k);
                if (bnd->s + kc->ext < end) {
                    //P; out of scope.
                    move_uniq(kc, e, bnd->s, end - 2);
                    Mantra copy = *bnd;
                    copy.e = end;
                    buf_grow_add(kc->bnd, 1ul, copy);
                    buf_grow_ks(kc, e, 1, &k);
                    excise_one(kc, e, k);
                } else {
                    //P; in scope of previous unique or start; excision
                    ++k;
                    buf_grow_ks(kc, e, (k - e->tgtk) - e->moved, &k);
                    for (uint32_t *q = e->tgtk + e->moved; q < k; ++q)
                        excise_one(kc, e, q);
                    --k;
                }
                bnd->s = end + 2;
            }
            ++k;
        }

        if (bnd->s + kc->ext >= bnd->e) {
            //P; in scope of previous unique or start; excision
            // (no bnd copy from e->obnd to kc->bnd)
            buf_grow_ks(kc, e, (k - e->tgtk) - e->moved, &k);
            for (uint32_t *q = e->tgtk + e->moved; q < k; ++q)
                excise_one(kc, e, q);
        } else {
            //P; alt
            move_uniq(kc, e, bnd->s, bnd->e - 2);
            buf_grow_add(kc->bnd, 1ul, *bnd);//K;
        }
    } while (++bnd != e->obnd + e->obnd_l);

    for (; h - kc->h != kc->h_l; ++h)
        hkoffs = update_append_hkoffs(kc, e, hkoffs);

    denser_k(kc, e, kc->hkoffs + kc->h_l, k);//K;
    //B; final state
}

// not valid for same address or with side effects
#define _TSWAP(t, a, b) do { t=a; a=b; b=t; } while(0)
#define _XSWAP(a, b) (a^=b, b^=a, a^=b)

static int
extd_uniqbnd(Key_t *kc)
{
    Iter_t e = {0};

    // should ideally not stop at a certain readlength, but rather when there are
    // no more dups, but then this goes on and on, for hg38 beyond extension 1845.

    unsigned end = (kc->readlength - KEY_WIDTH + 1) << 1;
    kc->ext_iter = buf_init(kc->ext_iter, 1);
    for (kc->ext = 0; kc->hkoffs[kc->h_l-1] && kc->ext != end; kc->ext += 2) {
        unsigned iter = 0;
        do { // until no no more new uniques
            Mantra *bnd = e.obnd;
            unsigned skctl = kc->kct_l;
            e.obnd = kc->bnd;
            kc->bnd = bnd;

            // The swapping is tricky. _l needs to be overwritten, _m needs to be
            // swapped to ensure the 1st allocation occurs correctly.
            e.obnd_l = kc->bnd_l;
            _XSWAP(e.obnd_m, kc->bnd_m);

            buf_grow(kc->bnd, 0); /* reserve as much space as in previous iteration */
            kc->bnd_l = 0;        /* ..but consider it as uninitialized: */

            ext_uq_iter(kc, &e);
            NB(skctl == kc->kct_l, "skctl(%u) != kc->kct_l(%u)", skctl, kc->kct_l);
            EPR("observed %u potential in iteration %u, extension %u\n",
                kc->ct, ++iter, kc->ext >> 1);
        } while (kc->ct > 0);

        EPR("----[ end of extension %u ]-------", kc->ext >> 1);
        buf_grow_add(kc->ext_iter, 1ul, iter);
    }
err:
    free(e.obnd);
    return 0;
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
            EPR("%s file present, refusing to overwrite.", ext[1+mode]);
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


