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
    char C* hdr = kc->id + h->part[0];

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
    ASSERT(p < (kc->s_l << 2), return -EFAULT, Pfmt "/%lu", p, kc->s_l);
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
        EPQ(IS_DBG_K(kc, k), "excised: [%lx]", K_OFFS(kc, k));
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

EPQ0(dbg >5, "[" Pfmt "%c]:\t", b2pos, IS_UQ(kct) ? '*' : ' ');
print_dna(seq.dna, dbg >5);

        _EVAL(get_nextnt(kc, p));
        seq.dna = _seq_next(res, seq);
        seq.p = _B2POS_OF(kc, kct);
        if (seq.p >= p++) { // first occurance (pos not yet reset)

            *kct ^= ((uint64_t)!seq.t << ORIENT_SHFT) ^ seq.p ^ p;

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
            EPQ(IS_DBG_K(kc, kct), "[" Pfmt "]", p);
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
        EPR("b2pos != b2end: " Pfmt " != " Pfmt " (happens for Y on hg19", b2pos, b2end + 1);
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
        EPR0("following pos " Pfmt " is: ", b2pos);
        for(i = 0, --p; i != KEY_WIDTH; ++i, ++p)
             EPR0("%c", b6(((kc->s[p>>2] >> ((p&3) << 1)) & 3) << 1));
        EPR0("\n[b2pos:" Pfmt ", koffs:%lx], \t", b2pos, K_OFFS(kc, kct));
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

// reset and remove non-uniqs, move uqs to start of kc->kct array
static inline int
extd_uq_by_k(kct_t* kc, pos_t p, pos_t pend, uint64_t *k, std::set<uint64_t> &pk)
{
    uint64_t infior;
    seq_t *ndxkct;
    int res;
    keyseq_t seq = {0};

    // of uniques the position is correct.
    EPQ(dbg > 4, Pfmt ".." Pfmt " excised", p, pend);

    ASSERT(pend < (kc->s_l << 2), return -EINVAL, Pfmt "/%lu?", pend, kc->s_l);
    ASSERT(p < pend, return _PRNT_SEQ_BY_POS(kc, pend), Pfmt " >= " Pfmt "?", p, pend);

    if (p + kc->ext + 1 >= pend && p + 1 != pend) { // within scope and keys to re-evaluate

        // could continue with key if adjacent
        seq.p = p;
        p -= KEY_WIDTH;
        _build_key(kc, seq, p, seq.p, seq.t);
        ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);

        // TODO: extension - KEY_WIDT as high bits of inferiority.
        // low inferiority bits zero if uniques were of previous extensions
        infior = kc->kct[*ndxkct] & MAX_INFIOR;
        if (k && (*k & MAX_INFIOR) > infior)
            infior = *k & MAX_INFIOR;

        expect(infior < MAX_INFIOR) {
            infior += INFIOR;
        } else {
            WARN("MAX_INFIOR reached");
        }

//        print_seq(&seq);
//        EPR("DECR, " Pfmt "\t" Pfmt, p, pend);
        do {
            _EVAL(get_nextnt(kc, p));
            seq.dna = _seq_next(res, seq);

            ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
            k = kc->kct + *ndxkct;
//            EPR("..%lu..", _B2POS_OF(kc, k));
//            print_seq(&seq);
            update_max_infior(k, infior);

            // in next iteration: reeval whether this is still a dup.
            //ASSERT(k > kc->kct_scope[0], return print_seq(&seq),
            //        "DUPBIT unset before kct occurance?");
            if (*k & DUP_BIT) {
                if (_B2POS_OF(kc, k) > p) {
                    ++kc->reeval;
                } else if (_B2POS_OF(kc, k) == p) {
                    // all remaining were excised this iteration
                    --kc->reeval; // this key won't extend
                    *k = 0; // position is not certain, but it happend in this iteration.
                }
            } else {
                ASSERT(_B2POS_OF(kc, k) < p, res = -EFAULT; goto err);
                // not certain yet, but may be unique. They are then handled in next iter.
            }

            if (IS_DBG_K(kc, k)) dbgndxkct = *ndxkct;
            EPQ(IS_DBG_K(kc, k), "excised: [%lx, %lu, " Sfmt "]", K_OFFS(kc, k), B2POS_OF(*k), *ndxkct);
        } while (++p != pend);

        ASSERT(IS_UQ(k), return -EFAULT, Pfmt ", %lu", p, _B2POS_OF(kc, k));

        // need to determine kcts for inbetween. Same can occur twice, use end instead.
    }
    res = 0;
err:
    return res;
}

static void
track_regions(kct_t* kc, std::list<Hdr*>::iterator& h, C pos_t p)
{
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

}

static int
reached_boundary(kct_t* kc, std::list<Hdr*>::iterator& h, C pos_t prev, C pos_t p, std::set<uint64_t> &pk)
{
    int res = 0;
    C pos_t pend = (*kc->bdit).e + (*h)->s_s;
    C pos_t b2start = (*kc->bdit).s + (*h)->s_s + KEY_WIDTH;
    if ((prev + kc->ext > pend) && (pend > prev)) {
        if ((*kc->bdit).s == b2start) {
            kc->bdit = (*h)->bnd.erase(kc->bdit); EPQ(dbg > 4, "(*h)->bnd.erase(kc->bdit)");
        } else {
            ++kc->bdit; EPQ(dbg > 4, "++kc->bdit(1)");
        }
        _EVAL(extd_uq_by_k(kc, prev, pend, NULL, pk));

    } else {
        (*kc->bdit).e = pend - KEY_WIDTH;
        ++kc->bdit; EPQ(dbg > 4, "++kc->bdit(2) => possibly too late: occurs only at uniq (only has a valid pos)");
        if ((dbg > 1) && ((p - prev) > 1000000)) {
            keyseq_t seq = {0};
            seq.p = p - KEY_WIDTH;
            _build_key(kc, seq, seq.p, p, seq.t);
            seq_t ndx;
            ndx = _get_kct0(kc, seq, seq.t, ndx, return -EFAULT);
            EPR("1M+ diatance between uniques " Sfmt, ndx);
            _PRNT_SEQ_BY_POS(kc, prev);
            _PRNT_SEQ_BY_POS(kc, p);
        }
    }
    track_regions(kc, h, p);
err:
    return res;    
}

static int emplace_kct(kct_t* kc,  uint64_t *sk,  uint64_t *k)
{
    keyseq_t seq = {0};
    seq.p = _B2POS_OF(kc, k);
    pos_t p = seq.p - KEY_WIDTH;
    _build_key(kc, seq, p, seq.p, seq.t);
    seq_t *ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);

    // ndxkct points to reordered kcts
    *ndxkct = k - kc->kct;
    
    // also swap info at kcts
    *sk = *k;
    return 0;
}

static int swap_kct(kct_t* kc,  uint64_t *k1,  uint64_t *k2, seq_t *ndxkct2, uint32_t line)
{
    EPQ(k1 == k2, "k1 == k2");

//    EPR("swap %lu <-> %lu => uq", _B2POS_OF(kc, k1), _B2POS_OF(kc, k2));
    // calc pos (to seq guarantee) => dna + rc => ndx and also swap ndxct!

    // swap ndxcts: TODO: one of these may not have to be recalculated in next iter.
    keyseq_t seq = {0};
    seq.p = _B2POS_OF(kc, k1);
    pos_t p = seq.p - KEY_WIDTH;
    _build_key(kc, seq, p, seq.p, seq.t);
    seq_t *ndxkct1 = _get_kct(kc, seq, seq.t, return -EFAULT);

    if (dbg > 4) {
        EPR("swap_kct() called on %s:%u", __FILE__, line);
        EPR("swap%u in %lu .. %lu:", ndxkct2 == NULL, B2POS_OF(*k1), B2POS_OF(*k2));
    }

    if (dbg > 4 || dbgndxkct == *ndxkct1 || dbgndxkct == *ndxkct2) {
        EPR("--");
        if (dbgndxkct == *ndxkct1 || dbgndxkct == *ndxkct2) {
            dbgndx = (dbgndxkct == *ndxkct1 ? ndxkct2 : ndxkct1) - kc->ndxkct;
            dbgk = *k2;
        }
        _PRNT_SEQ_BY_POS(kc, seq.p);
        _PRNT_SEQ_BY_POS(kc, B2POS_OF(*k2));
    }

    // ndxkct points to reordered kcts
    uint64_t t = *ndxkct1;
    *ndxkct1 = *ndxkct2;
    *ndxkct2 = t;
    
    // also swap info at kcts
    t = *k1;
    *k1 = *k2;
    *k2 = t;
    //ASSERT(p == _B2POS_OF(kc->kct + *ndxkct2), return -EFAULT);
    return 0;
}

static int
place_uniques(kct_t* kc, uint64_t* sk, uint64_t* &kend, std::set<uint64_t> &pk)
{
    pos_t p;
    int res;

    for (std::set<uint64_t>::iterator it = pk.begin(); it != pk.end(); ++it) {
        uint64_t kv = *it;
        if (~kv & DUP_BIT) {
            ++kc->uqct;
            p = B2POS_OF(kv);
            keyseq_t seq = {0};
            seq.p = p - KEY_WIDTH;
            _build_key(kc, seq, seq.p, p, seq.t);
            seq_t *ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
            // If smaller than sk, an ambiguous key was apparently unique. These are handled in
            // next iteration
            ASSERT(kc->kct + *ndxkct >= sk, return -EFAULT);

            _EVAL(swap_kct(kc, kc->kct + *ndxkct, kend--, ndxkct, __LINE__));
        }
    }

err:
    return res;
}

// FIXME/TODO/XXX: loop kc->ext from 1 to readlength
// order of kcts is important, possibly also for not yet uniq, maintain?

static int
ext_uq_iter(kct_t* kc)
{
    int res;

    uint64_t *kend = kc->kct + kc->kct_l - kc->last_uqct - 1; // location after uniques, from last time
    std::list<Hdr*>::iterator h = kc->h.begin();
    ASSERT(kc->uqct < kc->kct_l, return -EFAULT);
    std::set<uint64_t> pk;     // storage for unique key offsetss
    kc->bdit = (*h)->bnd.begin();
    EPR("[%s] %lu .. %lu", kc->id + (*h)->part[0], (*kc->bdit).s + (*h)->s_s,
            (*kc->bdit).e + (*h)->s_s);
    pos_t b2start = (*kc->bdit).s + KEY_WIDTH;
    pos_t b2end = (*kc->bdit).e;
    pos_t p, prev = b2start - 1;
    start_dbg(kc, *h, 0ul);
    // iterate over keys
    // 1) reset and remove non-uniqs, move uqs to start of array
    // 2) determine, extend unique regions and decr_excise *only new* ones
    //    new uqs lie within mantras
    // 3) shrink mantras accordingly.
    //
    // Only keys on mantra are processed (prev. processed lie before start k pointer).
    //

    // keys are kept sorted, first on uniqueness, then on position [TODO for non-uniques].

    //if (kc->iter > 0)
    //    dbg = 6;

    // dna ^ rc => ndx; ndxct[ndx] => kct => pos (regardless of whether uniq: s[pos] => dna and rc)
    uint64_t *sk = kc->kct;
    for (uint64_t *k = kc->kct; k < kend; ++k) {
        seq_t *ndxkct = NULL;

        if (dbg > 5)  _PRNT_SEQ_BY_POS(kc, B2POS_OF(*k));
        EPQ(IS_DBG_K(kc, k), "debug k: [%lx, %lu]", K_OFFS(kc, k), B2POS_OF(*k)); // 1
        if (IS_UQ(k)) {
            p = _B2POS_OF(kc, k);
            EPQ(dbg > 4, "%lu%s", B2POS_OF(*k), p - prev < kc->ext + 1 ? " .. 2nd uq" : " => 1st uq");
            pk.insert(*k);
            keyseq_t seq = {0};
            ASSERT(p != prev || prev == b2start, return -EFAULT, "%lu", pk.size());
            if (p > b2end + (*h)->s_s) {

                // handle last boundary
                _EVAL(reached_boundary(kc, h, prev, p, pk));
                //ASSERT(p >= (*kc->bdit).s + (*h)->s_s + KEY_WIDTH, return -EFAULT,
                //    Pfmt ", %lu", p, (*kc->bdit).s + (*h)->s_s + KEY_WIDTH);
                b2end = (*kc->bdit).e;
                b2start = (*kc->bdit).s;
                EPQ(kc->iter == 0 || dbg > 5, "[%s] %lu .. %lu",
                        kc->id + (*h)->part[0], b2start + (*h)->s_s, b2end + (*h)->s_s);
                b2start += KEY_WIDTH;
                prev = b2start - 1;
            }

            if (p - prev < kc->ext + 1) { // a 2nd uq
                if (prev == b2start - 1) // adjoining boundary
                    (*kc->bdit).e = p;
                else
                    ASSERT (prev < p, return -EFAULT);

                // what if all remaining occurances happen between two boundaries?
                if (prev - p > 1)
                    _EVAL(extd_uq_by_k(kc, prev, p, k, pk));
                prev = p;
                if (!ndxkct) {
                    seq.p = _B2POS_OF(kc, k);
                    p = seq.p - KEY_WIDTH;
                    _build_key(kc, seq, p, seq.p, seq.t);
                    ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
                } else {
                    p = _B2POS_OF(kc, k);
                    uint32_t* test = _get_kct(kc, seq, seq.t, return -EFAULT);
                    ASSERT(p == _B2POS_OF(kc, kc->kct + *ndxkct), return -EFAULT);
                    ASSERT(ndxkct == test, return -EFAULT);
                }
            } else { // first occurance. need to close previous, if set.


                //EPR("update: " Pfmt " .. " Pfmt, prev, p);
                // 1: from prev till now, add position and dup reevaluation
                seq.p = prev + 1 - KEY_WIDTH;
                _build_key(kc, seq, seq.p, prev + 1, seq.t);
                ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
                do {
                    k = kc->kct + *ndxkct;
                    EPQ(IS_DBG_K(kc, k), "debug k (2): [%lx, %lu, " Sfmt "]",
                            K_OFFS(kc, k), B2POS_OF(*k), *ndxkct);
                    // XXX: for this to work also non-uniq kcts should be sorted on pos?
                    if (*k & DUP_BIT) {
                        if (_B2POS_OF(kc, k) > seq.p) { // first occurance of pot.multiple mv to start.
                            *k &= ~DUP_BIT; // unset for 1st occurance
                            ASSERT(sk <= k, return -EFAULT);
                            _EVAL(swap_kct(kc, sk, k, ndxkct, __LINE__));
                            ++kc->reeval;
                            k = sk++;
                        } else if (_B2POS_OF(kc, k) == seq.p) {
                            // all but last positions were excised, unique after all!
                            pk.insert(*k);
                            unsigned len = seq.p - prev;
                            if (len < kc->ext + 1) { EPQ(dbg > 5, "rollback");
                                if (len > 1)
                                    _EVAL(extd_uq_by_k(kc, prev, seq.p, k, pk));
                                prev = seq.p;
                                if (p - prev < kc->ext + 1) { EPQ(dbg > 3, "rollback entirely");
                                    seq.p = p - KEY_WIDTH;
                                    _build_key(kc, seq, seq.p, p, seq.t);
                                    ndxkct = _get_kct(kc, seq, seq.t, return -EFAULT);
                                    if (p - prev > 1)
                                        _EVAL(extd_uq_by_k(kc, prev, p, kc->kct + *ndxkct, pk));
                                    prev = p;
                                    break;
                                }
                            }
                        }
                    } else if (_B2POS_OF(kc, k) < seq.p) { // no dup after all
                        *k |= DUP_BIT;
                        --kc->reeval;
                    }
                    *k &= INFIOR_MASK; // unset strand bit and pos (dupbit is highest and preserved)
                    *k |= ((uint64_t)!seq.t << ORIENT_SHFT) | seq.p; // set new pos and strand
                    _EVAL(get_nextnt(kc, seq.p));
                    seq.dna = _seq_next(res, seq);
                    ndxkct = _get_kct(kc, seq, seq.t, res = -EFAULT; goto err);
                } while (++seq.p != p);

                if (prev != p) {
                    // a new region starts if not adjoining end or at start
                    if (p <= b2end + (*h)->s_s && prev != b2start - 1)
                        start_region(kc, *h, b2start, prev);
                    // a region ends if not adjoining start or at end
                    if (b2start + (*h)->s_s <= p - KEY_WIDTH && p != b2end)
                        (*kc->bdit).e = p - KEY_WIDTH;
                    prev = p;
                } // else rollback entirely occurred.
                k = kc->kct + *ndxkct;
            }
        }
    }
    // last region & excision ahould already have been occurred.
    _EVAL(place_uniques(kc, sk, kend, pk));
    //ASSERT(0, return -EFAULT, "(NOT!;) end of loop reached! (TODO: swapping..)");
    // FIXME: iteration over swapped instead of uniques
    //

    // TODO: some uniques no longer occur: all are excised.
    kc->last_uqct = kc->uqct;

    /*for (h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
    }*/
    EPQ(dbg > 0, "observed %u uniques in iteration %u, extension %u, %u to be reevaluated\n",
            kc->uqct, ++kc->iter, kc->ext, kc->reeval);

    res = kc->reeval;
err:
    if (res < -1) {
        EPR("prev:" Pfmt ", p:" Pfmt ", b2end:" Pfmt, prev, p, b2end);
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


