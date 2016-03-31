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
#include <queue>
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"

static int
start_dbg(kct_t *C kc, Hdr *C h, seq_t C dna)
{
    char C* hdr = kc->id + h->part[0];

    // TODO: if not start of contig or just after a N-stretch, the dna key should be
    // just after an unique key. we could check this.

    dbg = strncmp(hdr, dbgchr, strlen(hdr)) ? 3 : 5;

    if (dbg > 3) {
        EPR0("----[\t%s%s:(%u+)%u..%u\t]----\tdna:", strlen(hdr) < 8 ? "\t": "",
                hdr, (*kc->bdit).corr, (*kc->bdit).s, (*kc->bdit).e);
        print_dna(dna);
        show_mantras(kc, h);
    }
    return 0;
}

void
free_kc(kct_t* kc)
{
    std::list<Hdr*>::iterator it;
    for (it = kc->h.begin(); it != kc->h.end(); ++it) {
        Hdr* h = *it;
        free(h->part);
        delete h;
    }
    _buf_free(kc->id);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->ndxkct);
}

static inline void
update_max_infior(uint64_t *C k, uint64_t infior)
{
    expect(infior < MAX_INFIOR) {
        infior += INFIOR;
    } else {
        WARN("MAX_INFIOR reached");
    }
    if (infior > *k)
        *k ^= (*k ^ infior) & INFIOR_MASK;
}

static inline C int
get_nextnt(kct_t C*C kc, pos_t p)
{
    ASSERT(p < (kc->s_l << 2), return -EFAULT, "%lu/%lu", p, kc->s_l);
    return (kc->s[p>>2] >> ((p&3) << 1)) & 3;
}

static void
decr_excise_all(kct_t *C kc, uint64_t* k, uint64_t C*C kct, unsigned const end)
{
    if (k == NULL) return;
    uint64_t infior = (kct ? std::max(*k, *kct) : *k) & MAX_INFIOR;
    expect(infior < MAX_INFIOR) {
        infior += INFIOR;
    } else {
        WARN("MAX_INFIOR reached");
    }
    unsigned i = IS_UQ(k) ? 1 : 0;
    while (i < end) {
        uint64_t *C k = kc->kct_scope[i++];
        if (k != kct) // .. don't elevate key if it is the one that became uniq
            update_max_infior(k, infior);
        EPQ(IS_DBG_K(kc, k), "exiced: [%u]", K_OFFS(kc, k));
        ++kc->uqct;
    }
}

static inline void start_region(kct_t* kc, Hdr *C h, C pos_t b2start, C pos_t last_uq_pos)
{
    if ((*kc->bdit).e == b2start) {
        (*kc->bdit).s = last_uq_pos - KEY_WIDTH;
    }  else if ((*kc->bdit).e != last_uq_pos) { // not downstream adjoining.
        h->bnd.insert(kc->bdit, *kc->bdit); // insert a copy of the current
        // start is earlier to enable building of the 1st key.
        (*kc->bdit).s = last_uq_pos - KEY_WIDTH;
    }
}

/*
 * A first uniq marks a potential start of a region, a 2nd or later uniq, within
 * scope - keylength minus keywidth distance - triggers a region insertion or
 * update. Also retrieve for this index the offset.
 */
        // TODO: if all nextNts are the same increase keylength.

// reverse seq could be required for mapping - maybe 
// N-stretches, later also skips, future: SNVs, splice sites?
static int
ext_uq_bnd(kct_t *C kc, Hdr *C h)
{
    uint64_t *kct = NULL, *k;
    unsigned i, index = 0, ext = kc->ext;
    C uint32_t b2end = (*kc->bdit).e;
    C uint32_t b2start = (*kc->bdit).s;
    keyseq_t seq = {0};
    uint32_t b2pos = b2start;
    int res;

    // build 1st key
    pos_t p = b2pos + h->s_s;
    seq.p = p + KEY_WIDTH;
    ASSERT_SCOPE_RNG(kc, p, seq.p, return -EFAULT);
    _build_key(kc, seq, p, seq.p, seq.t);
    start_dbg(kc, h, seq.dna);

    kc->kct_scope[0] = NULL;

    for (b2pos += KEY_WIDTH; b2pos <= b2end; ++b2pos) { // until next uniq region, stretch or contig
        /* process new key */
        kct = kc->kct + *(_get_kct(kc, seq, seq.t, res = -EFAULT; goto err));
        if (IS_DBG_K(kc, kct)) {
            DESCRIBE_KEY(kc, kct, '>');
            print_seq(&seq, IS_DBG_K(kc, kct));
        }

        ++h->total;
EPQ0(dbg >5, "[%u%c]:\t", b2pos, IS_UQ(kct) ? '*' : ' ');
print_dna(seq.dna, dbg >5);

        _EVAL(get_nextnt(kc, p));
        seq.dna = _seq_next(res, seq);
        seq.p = B2POS_OF(kct);
        if (seq.p >= p++) { // first occurance (pos not yet reset)

            *kct ^= (!seq.t << ORIENT_SHFT) ^ seq.p ^ p;

	    k = kc->kct_scope[0];
            if (k == NULL || (b2start + index + KEY_WIDTH - 1) == b2pos) {
                decr_excise_all(kc, k, kct, index);
                (*kc->bdit).e = b2start;

            } else if (IS_UQ(k)) {

                decr_excise_all(kc, k, kct, index);

            } else {
                (*kc->bdit).e = b2pos;
            }
	    index = 0;

        } else {

            *kct |= DUP_BIT;
            // TODO: count all same Nts and store seq, pos or ref to these.
            // strand bit (orientation first pos) is no longer functional:
            // can indicate whether seq & len is stored, or a position
            EPQ(IS_DBG_K(kc, kct), "[%u]", p);
            if (index == ext + 1) {
                i = 0;
                if (IS_UQ(kc->kct_scope[i])) {
                    start_region(kc, h, b2start, b2pos - index);
                    ++i;
                }
                index = 0;
            }
        }
        ASSERT(index <= ext, return -EFAULT);
        kc->kct_scope[index++] = kct;
    }
    EPQ(IS_DBG_K(kc, kct), "last key was dbgtsoffs");
    if (b2pos != b2end + 1) {
        EPR("b2pos != b2end: %u != %u (happens for Y on hg19", b2pos, b2end + 1);
        show_mantras(kc, h);
        ASSERT(strncmp(kc->id + h->part[0], "Y", strlen("Y")) == 0, return -EFAULT);
    }
    k = kc->kct_scope[0];
    ASSERT(k != NULL, res = -EFAULT; goto err);
    if (IS_UQ(k)) {
        if ((*kc->bdit).s == b2start) {
            kc->bdit = h->bnd.erase(kc->bdit);
        } else if (((*kc->bdit).e + index) == b2pos) {
            kc->bdit++;
        } else {
            kc->bdit++;
            //show_mantras(kc, h);
        }
        decr_excise_all(kc, k, NULL, index);

    } else {
        EPQ(dbg > 3, "(post_loop non-uq)");
//DESCRIBE_KEY(kc, k, 'k');
        (*kc->bdit).e = b2pos;
        ++kc->bdit;
    }
//show_mantras(kc, h);
    res = 0;
err:
    if (res < -1) {
        EPR0("following pos %u is: ", b2pos);
        for(i = 0, --p; i != KEY_WIDTH; ++i, ++p)
             EPR0("%c", b6(((kc->s[p>>2] >> ((p&3) << 1)) & 3) << 1));
        EPR0("\n[b2pos:%u, koffs:%lu], \t", b2pos, K_OFFS(kc, kct));
        print_seq(&seq);
        show_mantras(kc, h);
    }
    return res;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    //FIXME: rather than % of all display % mapable of remaining (mantras)
    kc->bdit = h->bnd.begin();
    IFOUT(kc->bdit == h->bnd.end(), "Already fully mapable: %s", kc->id + h->part[0]);

    EPQ(dbg > 3, "Processing %s", kc->id + h->part[0]);

    do {
        int ret = ext_uq_bnd(kc, h);
        if (ret < 0) return ret;
    } while (kc->bdit != h->bnd.end());

    IFOUT(h->bnd.begin() == h->bnd.end(), "Became fully mapable: %s", kc->id + h->part[0]);
#ifdef DEBUG
    if (dbg > 4)
        show_mantras(kc, h);
#endif
    EPQ(dbg > 2, "Contig %s", kc->id + h->part[0]);
out:
    return 0;
}

static inline int print_seq_from_pos(kct_t* kc, pos_t p)
{
    keyseq_t seq = {0};
    EPR0("%lu\t", p);
    seq.p = p;
    p -= KEY_WIDTH;
    _build_key(kc, seq, p, seq.p, seq.t);
    return print_seq(&seq);
}

// reset and remove non-uniqs, move uqs to start of kc->kct array
// FIXME: mantras?
static inline int
extd_uq_by_k(kct_t* kc, pos_t p, pos_t pend, uint64_t *k)
{
    uint64_t infior, t;
    uint32_t *ndxkct;
    int res;
    keyseq_t seq = {0};

    // of uniques the position is correct.
    EPQ(dbg > 4, "%lu..%lu exised", p, pend);

    ASSERT(pend < (kc->s_l << 2), return -EINVAL, "%lu/%lu?", pend, kc->s_l);
    ASSERT(p < pend, return print_seq_from_pos(kc, pend), "%lu >= %lu?", p, pend);

    if (p + kc->ext + 1 >= pend && p + 1 != pend) { // within scope and keys to re-evaluate

        // could continue with key if adjacent
        seq.p = p;
        p -= KEY_WIDTH;
        _build_key(kc, seq, p, seq.p, seq.t);
        ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);

        infior = kc->kct[*ndxkct] & MAX_INFIOR;
        if (k && (*k & MAX_INFIOR) > infior)
            infior = *k & MAX_INFIOR;

        expect(infior < MAX_INFIOR) {
            infior += INFIOR;
        } else {
            WARN("MAX_INFIOR reached");
        }

//        print_seq(&seq);
//        EPR("DECR, %lu\t%lu", p, pend);
        do {
            _EVAL(get_nextnt(kc, p));
            seq.dna = _seq_next(res, seq);

            ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
            k = kc->kct + *ndxkct;
//            EPR("..%lu..", B2POS_OF(k));
//            print_seq(&seq);
            update_max_infior(k, infior);

            // in next iteration: reeval whether this is still a dup.
            //ASSERT(k > kc->kct_scope[0], return print_seq(&seq),
            //        "DUPBIT unset before kct occurance?");
            if ((*k & DUP_BIT) && B2POS_OF(k) >= p) {
                *k &= ~DUP_BIT;
                ++kc->reeval;
            }

            EPQ(IS_DBG_K(kc, k), "exiced: [%u]", K_OFFS(kc, k));
        } while (++p != pend);

        ASSERT(IS_UQ(k), return -EFAULT, "%lu, %lu", p, B2POS_OF(k));
        // need to determine kcts for inbetween. Same can occur twice, use end instead.
    }
    res = 0;
err:
    return res;
}

static int
reached_boundary(kct_t* kc, std::list<Hdr*>::iterator& h, C pos_t lastp, C pos_t p)
{
    int res = 0;
    C pos_t pend = (*kc->bdit).e + (*h)->s_s;
    C pos_t b2start = (*kc->bdit).s + (*h)->s_s + KEY_WIDTH;
    if ((lastp + kc->ext > pend) && (pend > lastp)) {
        if ((*kc->bdit).s == b2start) {
            kc->bdit = (*h)->bnd.erase(kc->bdit); EPQ(dbg > 4, "(*h)->bnd.erase(kc->bdit)");
        } else {
            ++kc->bdit; EPQ(dbg > 4, "++kc->bdit(1)");
        }
        _EVAL(extd_uq_by_k(kc, lastp, pend, NULL));

    } else {
        (*kc->bdit).e = pend;
        ++kc->bdit; EPQ(dbg > 4, "++kc->bdit(2)");
    }
    for (;;) {
        if (kc->bdit == (*h)->bnd.end()) {
            if (dbg > 3)
                show_mantras(kc, *h);
            kc->bdit = (*++h)->bnd.begin(); EPQ(dbg > 4, "kc->bdit = (*++h)->bnd.begin()");
            start_dbg(kc, *h, 0ul);
        }
        if (p <= (*kc->bdit).e + (*h)->s_s)
            break;
        ++kc->bdit;  EPQ(dbg > 4, "++kc->bdit(3)");
    }
err:
    return res;    
}
// keep track of the mantra we're updating currently. With this, every processed key should
// lie between (*kc->bdit).s and (*kc->bdit).e. Enables extending the unique regions (and
// decreasing or even removing mantras accordingly).
static inline void
track_mantra(kct_t* kc, std::list<Hdr*>::iterator *h, C pos_t p)
{
    while (p > ((*kc->bdit).e + (**h)->s_s)) {
        EPR("bnd %u", p);
        if (++kc->bdit == (**h)->bnd.end()) {
            kc->bdit = (*++*h)->bnd.begin();
            EPR("%s", kc->id + (**h)->part[0]);
        }
    }
}

static int swap_kct(kct_t* kc,  uint64_t *uk,  uint64_t *nuk, uint32_t *ndxkct2)
{
    if (nuk == uk)
        return 0;
//    EPR("swap %lu <-> %lu => uq", B2POS_OF(uk), B2POS_OF(nuk));
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {0};
    seq.p = B2POS_OF(uk);
    pos_t p = seq.p - KEY_WIDTH;
    _build_key(kc, seq, p, seq.p, seq.t);
    uint32_t *ndxkct1 = _get_kct(kc, seq, seq.t, return -EFAULT);

    if (ndxkct2 == NULL) {
        seq.p = B2POS_OF(nuk);
        p = seq.p - KEY_WIDTH;
        _build_key(kc, seq, p, seq.p, seq.t);
        ndxkct2 = _get_kct(kc, seq, seq.t, return -EFAULT);
    }

    // ndxkct points to reordered kcts
    uint64_t t = *ndxkct1;
    *ndxkct1 = *ndxkct2;
    *ndxkct2 = t;
    
    // also swap info at kcts
    t = *uk;
    *uk = *nuk;
    *nuk = t;
    return 0;
}

// FIXME/TODO/XXX: loop kc->ext from 1 to readlength
// order of kcts is important, possibly also for not yet uniq, maintain?

static int
ext_uq_iter(kct_t* kc)
{
    uint64_t totNts = 0ul; // FIXME: put in kc and move to key_init
    int res;
    EPQ(dbg > 5, "Clearing dups for next iteration");

    uint64_t *nuk = kc->kct;            // location for unique keys
    uint64_t *sk = nuk;
    uint64_t *uk = nuk + kc->kct_l - kc->last_uqct - 1; // location after uniques, from last time
    std::list<Hdr*>::iterator h = kc->h.begin();
    ASSERT(kc->uqct < kc->kct_l, return -EFAULT);
    std::queue<uint64_t> pk;     // storage for unique/non-unique keys
    kc->bdit = (*h)->bnd.begin();
    EPR("[%s] %lu .. %lu", kc->id + (*h)->part[0], (*kc->bdit).s + (*h)->s_s, (*kc->bdit).e + (*h)->s_s);
    pos_t b2start = (*kc->bdit).s + KEY_WIDTH;
    pos_t b2end = (*kc->bdit).e;
    pos_t p, lastp = b2start;
    start_dbg(kc, *h, 0ul);
    // iterate over keys
    // 1) reset and remove non-uniqs, move uqs to start of array
    // 2) determine, extend unique regions and decr_excise *only new* ones
    //    new uqs lie within mantras
    // 3) shrink mantras accordingly.
    //
    // Only keys on mantra are processed (prev. processed lie before start nuk pointer).
    //

    // keys are kept sorted, first on uniqueness, then on position [TODO for non-uniques].


    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    while (nuk != uk) {
        ASSERT(uk - kc->kct < kc->kct_l && uk > 0, return -EFAULT, "uk");
        uint32_t *ndxkct = NULL;
        if (IS_UQ(nuk)) {
            ++kc->uqct;
            p = B2POS_OF(nuk);
            EPQ(dbg > 4, "%lu => uq", p); //
            if (p > b2end + (*h)->s_s) {

                // handle last boundary
                _EVAL(reached_boundary(kc, h, lastp, p));
                //ASSERT(p >= (*kc->bdit).s + (*h)->s_s + KEY_WIDTH, return -EFAULT,
                //    "%lu, %lu", p, (*kc->bdit).s + (*h)->s_s + KEY_WIDTH);
                b2end = (*kc->bdit).e;
                b2start = (*kc->bdit).s + (*h)->s_s;
                EPR("(%lu uq)", kc->uqct);
                EPR("[%s] %lu .. %lu", kc->id + (*h)->part[0],
                        b2start, b2end + (*h)->s_s);
                lastp = b2start += KEY_WIDTH;
            }

            if (lastp + kc->ext + 1 >= p) { // a 2nd uq
                if (lastp == b2start) // adjoining boundary
                    (*kc->bdit).e = b2start;

                if (lastp < p)
                    _EVAL(extd_uq_by_k(kc, lastp, p, nuk));
                lastp = p;
            } else { // first occurance. need to close previous, if set.
                // 1: mark start of a new region
                if ((*kc->bdit).e <= b2end)
                    start_region(kc, *h, b2start, lastp);
                (*kc->bdit).e = p;

                //EPR("update: %lu .. %lu", lastp, p);
                // 2: from lastp till now, add position and dup reevaluation
                keyseq_t seq = {0};
                seq.p = lastp + 1;
                lastp = p;
                p = seq.p - KEY_WIDTH;
                ASSERT(p < lastp, return -EFAULT);
                _build_key(kc, seq, p, seq.p, seq.t);
                ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);

                do {
                    uint64_t *k = kc->kct + *ndxkct;
                    // XXX: for this to work also non-uniq kcts should be sorted on pos?
                    if (B2POS_OF(k) >= p) {
                        if (*k & DUP_BIT) { // first occurance, move key to start
                            *k &= ~DUP_BIT; // don't set it yet at 1st occurance
                            _EVAL(swap_kct(kc, sk++, k, ndxkct));
                            ++kc->reeval;
                        } else {
                            *k |= DUP_BIT;
                            --kc->reeval; // no dup after all
                        }
                    }
                    *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
                    *k |= (!seq.t << ORIENT_SHFT) | p; // set new pos and strand
                    _EVAL(get_nextnt(kc, p));
                    seq.dna = _seq_next(res, seq);
                    ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
                } while (++p != lastp);

            }
            _EVAL(swap_kct(kc, uk--, nuk, ndxkct));
        }
        ++nuk;
    }

    kc->last_uqct = kc->uqct;

    /*for (h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
        totNts += (*h)->total;
    }*/

    EPQ(dbg > 0, "extended %u unique ranges in iteration %u, extension %u, %u to be reevaluated\n",
            kc->uqct, ++kc->iter, kc->ext, kc->reeval);

    res = kc->reeval;
err:
    if (res < -1) {
        EPR("lastp:%lu, p:%lu, b2end:%lu", lastp, p, b2end);
    }
    return res;
}

static int
extd_uniqbnd(kct_t* kc, struct gzfh_t* fhout)
{
    int res;
    kc->uqct = kc->reeval = 0;

    kc->kct_scope = (uint64_t**)malloc((kc->readlength - KEY_WIDTH+1) * sizeof(uint64_t*));
    ASSERT(kc->kct_scope != NULL, res = -EFAULT; goto err);
    
    for (kc->ext = res = 1; res > 0 && kc->ext != kc->readlength - KEY_WIDTH + 1;++kc->ext) {
        kc->iter = 0;
        do { // until no no more new uniques
            res = ext_uq_iter(kc);
        } while (res > 0);
    }
    _buf_free(kc->kct_scope);
    if (res == 0) {
        _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
        _ACTION(save_kc(fhout + 3, kc), "writing unique keycounts file");
    }
err:
    return res;
}

int
fa_index(struct gzfh_t* fh, unsigned readlength)
{
    int len, res = -ENOMEM;
    char file[1024];
    kct_t kc = {0};
    kc.readlength = readlength;
    unsigned mode;

    C char* ext[7] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub", ".uq"};

    if (fh[0].name) {
        len = strlen(fh[0].name);
    } else {
        ASSERT(fh[2].name != NULL, return -EFAULT);
        fh[0].name = file;

        len = strlen(fh[2].name);
        strncpy(fh[0].name, fh[2].name, ++len);
        _ACTION0(reopen(fh, ext[0], ext[5]), "")
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
         _ACTION0(reopen(fh, ext[5], ext[4]), "trying to open %s instead", ext[4])
        mode = fh[0].fp != NULL;
    }
    if (mode) {
        _ACTION(load_boundaries(fh, &kc), "loading boundary file %s", fh[0].name)
    }

    // keycount file did not exist, check twobit file.
    _ACTION0(reopen(fh, ext[mode == 2 ? 5 : 4], ext[1]),
                "%s does%s exist", ext[1], fh[0].fp ? "" : " not")
    _ACTION0(reopen(fh + 1, ext[5], ext[2]), "%s does%s exist", ext[2], fh[1].fp ? "" : " not")
    _ACTION0(reopen(fh + 3, ext[5], ext[3]), "%s does%s exist", ext[3], fh[3].fp ? "" : " not")

    if (mode && fh[0].fp && fh[1].fp && fh[3].fp) {
         _ACTION(load_seqb2(fh, &kc), "loading twobit sequence file")
         _ACTION(load_kc(fh + 3, &kc), "loading keycounts file")
    } else {
        bool found = false;
        for (int i=0; i != 4; ++i) {
            if (i == 2) continue;
            if (fh[i].fp) {
                found = true;
                EPR("%s file present, refusing to overwrite.", fh[i].name);
                rclose(fh +i);
            }
        }
        if (mode) {
            EPR("%s file present, refusing to overwrite.", ext[3+mode]);
            found = true;
        }
        if (found) goto err;
        EPR("starting from scratch.");
        _ACTION(fa_read(fh, &kc), "reading fasta")
    }
    if (mode < 2) {
        _ACTION(reopen(fh, ext[1], ext[5]), "")
        _ACTION(reopen(fh + 3, ext[3], ext[6]), "")
        _ACTION(extd_uniqbnd(&kc, fh), "extending unique boundaries")
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


