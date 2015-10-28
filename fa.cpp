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
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"

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
    _buf_free(kc->ts);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->ndxkct);
}

static inline void
update_max_infior(uint64_t *C kct, uint64_t infior)
{
    expect(infior < MAX_INFIOR) {
        infior += INFIOR;
    } else {
        WARN("MAX_INFIOR reached");
    }
    if (infior > *kct)
        *kct ^= (*kct ^ infior) & INFIOR_MASK;
}

static inline C int
check_key(kct_t C*C kc, C uint64_t cdna, uint64_t p)
{
    if (kc->s == NULL)
        return 0;

    uint64_t dna, rc, b, end = p + KEY_WIDTH;
    dna = rc = 0ul;

    for (; p != end; ++p) {
        dna = __rdsq(next, b, kc->s, p, dna, rc);
        EPR("%lu:%c", p, b6(b << 1));
    }

    if (dna == cdna) return 0;

    WARN("nextnt jump dna didn't match b2pos dna.");
    print_dna(cdna);
    print_2dna(dna, rc);

    return -EFAULT;
}

static inline C int
check_nt(kct_t C*C kc, C uint64_t nt, uint64_t p)
{
    if (kc->s == NULL)
        return 0;

    ASSERT(p < (kc->s_l << 2), return -EFAULT, "%lu/%lu", p, kc->s_l);
    p = (kc->s[p>>2] >> ((p&3) << 1)) & 3;
    expect(nt == p)
        return 0;

    WARN("nextnt didn't match twobit: sb2:%c, got %c", b6(p<<1), b6(nt<<1));
    return -EFAULT;
}

static inline C int
get_nextnt(kct_t C*C kc, uint64_t nt)
{
    nt &= B2POS_MASK;
    ASSERT(nt < (kc->ts_l << 2), return -EFAULT);
    return _GET_NEXT_NT(kc, nt); // next Nt
}

static inline int
mark_leaving(kct_t C*C kc, uint64_t *C k)
{
    int res = -EFAULT;
    uint64_t t, nt;
    EPQ(TSO_NT(k) == dbgtsoffs, "<%s%s%s%s%sREMAIN:%lu(-1?)\tNEXT_NT_NR:%lu\tPENDING:%lu",
            IS_FIRST(k) ? "FIRST\t" : "", ALL_SAME_NTS(k) ? "ALL_SAME\t" : "",
            IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",
            IS_LAST(k) ? "LAST\t" : "",REMAIN(k), FIRST_IN_SCOPE(k), PENDING(k));
    if (PENDING(k)) {
        *k -= ONE_PENDING;
        t = TSO_NT(k) + FIRST_IN_SCOPE(k); // position of nextnt, to compare with former nextnt.
        ASSERT(t < (kc->ts_l << 2), goto err);

        nt = _GET_NEXT_NT(kc, t);
        EPQ(TSO_NT(k) == dbgtsoffs, "\t%c", b6(nt<<1));

        if (IS_FIRST(k)) {
            if (IS_DISTINCT(k)) { // XXX: never taken
                EPQ(TSO_NT(k) == dbgtsoffs, " ==> distinct cleared");
                k[1] &= ~DISTINCT; // clear from last iteration
            } else if (kc->iter) { // XXX: never taken
                //if ok, then we can remove the kct iteration post loop?
//                ASSERT(k[1] & MARKED, goto err, "MARKED not obsolete"); // TODO: check after setting marked
            }
        } else {
            --t;
            EPQ(TSO_NT(k) == dbgtsoffs, "%lu:%c <=> %c", t, b6(_GET_NEXT_NT(kc, t)<<1), b6(nt<<1)); 
            ASSERT(t < (kc->ts_l << 2), goto err);
            if (IS_DISTINCT(k) == false && nt != _GET_NEXT_NT(kc, t)) {
                EPQ(TSO_NT(k) == dbgtsoffs, " ==> dbgtsoffs got distinct");
                k[1] |= DISTINCT;
            }
        }
        ++*k;
    }

    MARK_LAST(k)
    res = 0;
err:
    if (res || TSO_NT(k) == dbgtsoffs || dbgtsoffs == -2ul) {
        EPQ0(res, "[%lu]\t", TSO_NT(k));
        EPR("<%s%s%s%s%sREMAIN:%lu(-1?)\tNEXT_NT_NR:%lu\tPENDING:%lu",
                IS_FIRST(k) ? "FIRST\t" : "", ALL_SAME_NTS(k) ? "ALL_SAME\t" : "",
                IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",
                IS_LAST(k) ? "LAST\t" : "",REMAIN(k), FIRST_IN_SCOPE(k), PENDING(k));
    }
    return res;
}

static int
_mark_all(kct_t C*C kc, char C* f, unsigned l, unsigned i, unsigned const end)
{
    int res = 0;
    while (i < end && res == 0) {
        uint64_t* k = kc->kct_scope[i++];
        if (TSO_NT(k) == dbgtsoffs)
            EPR("%s +%u mark_leaving()", f, l);
        res = mark_leaving(kc, k);
    }
    if (res < 0)
        EPQ(res, "error during mark_all():%s +%u", f, l);
    return res;
}
#define mark_all(kc, i, e) _mark_all(kc, __FILE__, __LINE__, i, e)

static int
decr_excise(kct_t *C kc, uint64_t *C k)
{
    uint64_t nt, offs;
    uint8_t *q, *qe, c, t;
    int err = -EFAULT;

    if (ALL_SAME_NTS(k) == false) {// only when previously incremented
        ASSERT(PENDING(k), return -EFAULT);
        *k -= ONE_PENDING;
    }

    k[1] -= ONE_CT; // one less remains.
    ++kc->uqct;
    if (IS_LAST(k)) {//  src of movement would be target.
        err = 0;
        goto out;
    }

    offs = REMAIN(k), nt = TSO_NT(k);

    offs += nt; // remaining + ts offset => target of nt movement

    nt += LAST_IN_SCOPE(k); // add current position => src Nt (which is moved)

    q = &kc->ts[nt >> 2];
    qe = &kc->ts[offs >> 2];

    // can't mark k here since multiple same nt-keys may occur within scope.

    ASSERT((offs >> 2) < kc->ts_l, goto out, "nt:%lu\toffs*:%lu", nt, offs);
    ASSERT((nt >> 2) + 1 < kc->ts_l, goto out, "nt*:%lu\toffs:%lu", nt, offs);

    // excise out twobits, i.e. shift 2bits above current 2bit down.
    ASSERT(q <= qe, goto out, "nt:%lu > offs:%lu?", nt, offs);

    // shift for past Nt;
    c = ((nt & 3) + 1) << 1;
    t = *q & ((1ul << c) - 1ul); // capture this and past Nts.
    if (c) c -= 2;
    *q = ((*q ^ t) >> 2) | (*q & ((1ul << c) - 1ul)); // excise out Nt.
    t >>= c;
    EPQ(TSO_NT(k) == dbgtsoffs, "moving %c", b6(t << 1));

    // mask to cover this and past nucleotide.
    while (q != qe) {
        *q |= (q[1] & 3) << 6;
        *++q >>= 2;
    }
    // append excised Nt to end.
    offs = (offs & 3) << 1;        // set shift
    t <<= offs;                    // move excised in position
    offs = *q & ((1u << offs) - 1u); // below 2bits were shifted correctly
    *q = ((*q ^ offs) << 2) ^ t ^ offs;  // move top part back up, add rest.
    err = 0;
out:
    if (err || TSO_NT(k) == dbgtsoffs || dbgtsoffs == -2ul) {
        EPR("-%s%s%s%sREMAIN:%lu(-1?)\tNEXT_NT_NR:%lu\tPENDING:%lu",
                ALL_SAME_NTS(k) ? "ALL_SAME\t" : "",
                IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",
                IS_LAST(k) ? "LAST\t" : "",REMAIN(k), FIRST_IN_SCOPE(k), PENDING(k));
        EPQ(err, "[%lu]", TSO_NT(k));
    }
    return err;

}

static int
decr_excise_all(kct_t *C kc, uint64_t * k, uint64_t C*C kct, unsigned const end)
{
    int res = 0;
    uint64_t infior = kct ? max(*k, *kct) : *k;
    unsigned i = IS_UQ(k) ? 1 : 0; // XXX: never first taken.
    while (i < end) { // != test won't work for 1st if dummy.
        k = kc->kct_scope[i++];
        if (k != kct) // .. don't elevate key if it is the one that became uniq
            update_max_infior(k, infior);
        _EVAL(decr_excise(kc, k));
        MARK_LAST(k)
    }
err:
    EPQ(res, "%u/%u:before uniq %lu", i, end, TSO_NT(kct));
    return res;
}
#define decr_excise_all(kc, k, kct, end) ({\
    ASSERT((res = decr_excise_all(kc, k, kct, end)) >= 0, goto err,\
	"%s +%u decr_excise_all()", __FILE__, __LINE__);\
})

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

    C char* hdr = kc->id + h->part[0];

    uint64_t p, dna = (*kc->bdit).dna;   // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t *kct = NULL, *k;
    unsigned i, index = 0, ext = kc->ext;
    uint32_t b2pos;
    const uint32_t b2end = (*kc->bdit).e;
    int res;
    kc->kct_scope[0] = NULL;
    // Cannot check post jump dna, it is not stored in kc->s.

    dbg = strncmp(hdr, dbgchr, strlen(hdr)) ? 3 : 5;
    EPQ0(dbg > 3, "----[\t%s%s:(%u+)%u..%u\t]----\tdna:", strlen(hdr) < 8 ? "\t": "",
            hdr, (*kc->bdit).corr, (*kc->bdit).s, (*kc->bdit).e);
    print_dna(dna, dbg > 3);
    if (dbg > 3) show_mantras(kc, h);

    for (b2pos = (*kc->bdit).s; b2pos < b2end; ++b2pos) { // until next uniq region, stretch or contig
        uint64_t nt;
        p = b2pos + h->s_s;
        /* process new key */
        _get_new_kct(kc, kct, dna, rc);
        EPQ(TSO_NT(kct) == dbgtsoffs, "\t(%lu)", b2pos); // XXX

        if (IS_UQ(kct)) {

            *kct ^= (*kct ^ (p + 1)) & B2POS_MASK;
            _EVAL(get_nextnt(kc, kct[1]));
            nt = res;

            k = kc->kct_scope[0];
            if (k == NULL) {
                //everything already in place
            } else if (DOWNSTREAM_ADJOINING(index, b2pos)) {
                decr_excise_all(kc, k, kct, index);
                // mappable region extends beyond end of unique covered region, not before start.
                h->mapable += b2pos + ext - (*kc->bdit).s;
                (*kc->bdit).s = b2pos;
                (*kc->bdit).dna = dna;
            } else if (IS_UQ(k)) {
                if (FORMER_UQ_WAS_1ST(kc, index, b2pos)) {
                    _EVAL(insert_mantra(kc, h));
                    // the mantra ended and a new uniq region began
                    h->mapable -= (*kc->bdit).e - ext;
                } else {
                    h->mapable -= (*kc->bdit).s + ext;
                }

                decr_excise_all(kc, k, kct, index);
                //EPR("[%u, %u]", b2pos, index);
                h->mapable += b2pos + ext; // mapable region starts before 1st uniq.
                (*kc->bdit).s = b2pos;
                (*kc->bdit).dna = dna;

            } else { EPQ(dbg > 3, "first uniq");
                pot_mantra_end(kc, h, dna, b2pos);
                _EVAL(mark_all(kc, 0, index));
            }

	    index = 0;
        } else {

            if (ALL_SAME_NTS(kct)) {

                nt = _GET_NEXT_NT(kc, kct[1]);
                EPQ(TSO_NT(kct) == dbgtsoffs, "got %ul", nt);

            } else {
                _EVAL(get_nextnt(kc, NT_OFFS(kct)));
                nt = res;
                // IN_SCOPE (Nt's) must be instantly updated or next_nt won't match
                // twobit Nts for the same keys within scope.
                *kct += ONE_PENDING;
            }

	    if (index == ext) {
		i = 0;
		if (IS_UQ(kc->kct_scope[i])) {
                    if (FORMER_UQ_WAS_1ST(kc, index, b2pos)) { // XXX XXX XXX won't work
                        pot_mantra_end(kc, h, 0, 0); // undo (postpone) mantra ending
                    }
		    ++i;
                }

		_EVAL(mark_all(kc, i, index));
		index = 0;
            }
        }
        ASSERT(0 == check_nt(kc, nt, p & B2POS_MASK), res = -EFAULT; goto err, "%lu/%lu",
                TSO_NT(kct), NT_OFFS(kct));

        ASSERT(nt != -1ul, res = -EFAULT; goto err, "ndxkct 0x%lx [%u]", *kct, b2pos);
        ASSERT(index < ext, return -EFAULT);
        kc->kct_scope[index++] = kct;


        dna = _seq_next(nt, dna, rc);
    }
    EPQ(TSO_NT(kct) == dbgtsoffs, "last key was dbgtsoffs");
    ASSERT(b2pos == b2end, show_mantras(kc, h); res = -EFAULT; goto err, "[%u] %u", b2pos, b2end);
    k = kc->kct_scope[0];
    ASSERT(k != NULL, res = -EFAULT; goto err);
    if (IS_UQ(k)) {
        EPQ (dbg > 3, "Post loop adjoining boundary [%u,%u]", b2pos, h->mapable);
        if (FORMER_UQ_WAS_1ST(kc, index, b2pos)) {
            h->mapable -= (*kc->bdit).e - ext;
            kc->bdit++;
        } else {
            h->mapable -= (*kc->bdit).s + ext;
            kc->bdit = h->bnd.erase(kc->bdit);
        }
        decr_excise_all(kc, k, NULL, index);
        h->mapable += b2pos; // no extension beyond end

    } else if (k != NULL) {
	_EVAL(mark_all(kc, 0, index));
        EPQ(dbg > 3, "(post_loop3)");
        pot_mantra_end(kc, h, dna, b2pos);
        ++kc->bdit;
        if (kc->bdit != h->bnd.end() && (*kc->bdit).e <= b2pos) {
            //this seems necessary in case f N-stretch at end, why?
            kc->bdit = h->bnd.erase(kc->bdit);
        }
    }
    res = 0;
err:
    if (res < -1 || dbgtsoffs == -2ul) {
        EPR0("following is: ");
        for(i = 0; i != KEY_WIDTH; ++i, ++p)
             EPR0("%c", b6(((kc->s[p>>2] >> ((p&3) << 1)) & 3) << 1));
        EPR("");
        EPR0("[b2pos:%u, tsoffs:%lu], \t", b2pos, kct ? kct[1] & B2POS_MASK : ~0u);
        print_2dna(dna, rc);
    }
    return res;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    unsigned hdr_nts = h->end_pos + h->end_corr;
    kc->bdit = h->bnd.begin();
    IFOUT(kc->bdit == h->bnd.end(), "Already fully mapable: %s", kc->id + h->part[0]);

    EPR("Processing %s", kc->id + h->part[0]);

    h->mapable = 0ul;
    do {
        int ret = ext_uq_bnd(kc, h);
        if (ret < 0) return ret;
    } while (kc->bdit != h->bnd.end());

    IFOUT(h->bnd.begin() == h->bnd.end(), "Became fully mapable: %s", kc->id + h->part[0]);
#ifdef DEBUG
    if (dbg > 4)
        show_mantras(kc, h);
#endif
    EPQ(dbg > 2, "%s: %u/%u => %.2f%% mapable",
            kc->id + h->part[0], h->mapable, hdr_nts,
            h->end_pos ? 100.0f * h->mapable / hdr_nts : nanf("NAN"));
    ASSERT(h->mapable <= hdr_nts, show_mantras(kc, h); return -EFAULT);
out:
    return hdr_nts;
}

static int
ext_uq_iter(kct_t* kc)
{
    kc->uqct = kc->pending = 0u;
    uint64_t mapable = 0ul;
    uint64_t totNts = 0ul; // FIXME: put in kc and move to key_init
    for (std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
        mapable += (*h)->mapable;
        totNts += ret;
    }
    EPQ(dbg > 0, "extended %u unique ranges in iteration %u\n"
            "\t%lu/%lu => %.2f%% mapable. (%u pending)", kc->uqct, ++kc->iter,
            mapable, totNts, totNts ? 100.0f * mapable / totNts : nanf("NAN"), kc->pending);
    //dbg = 7;
    // FIXME: reset upon last occurance for non-uniq in inner loop
    for (uint64_t *k = kc->kct; k != kc->kct + kc->kct_l; k+=2) {
        if (IS_UQ(k) == false) {// at least one left (if unique position genomic 2bit)
#if defined(TEST_CLEARANCE)
            ASSERT(LAST_IN_SCOPE(k) == 0ul, return -EFAULT, "[%lu]", TSO_NT(k));
            //ASSERT(IS_DISTINCT(k) == false, return -EFAULT, "[%lu]", TSO_NT(k));
#else
            *k &= ~B2POS_MASK; // reset position tracking
            k[1] &= ~DISTINCT;
#endif
        }
    }

    return kc->uqct | kc->pending;
}

static int
extd_uniqbnd(kct_t* kc, struct gzfh_t** fhout)
{
    int res = -ENOMEM;
    size_t t;
    kc->iter = 0;

    // was ndx for storage, unset for pos, strand & infiority;
    for (unsigned i=0u; i != kc->kct_l; i += 2)
        kc->kct[i] = 0ul;

    t = kc->ext;
    kc->kct_scope = (uint64_t**)malloc(t * sizeof(uint64_t*));
    ASSERT(kc->kct_scope != NULL, res = -EFAULT; goto err);

    do { // until no no more new uniques
        res = ext_uq_iter(kc);
    } while (res > 0);
    _buf_free(kc->kct_scope);
    if (res == 0) {
        _ACTION(save_boundaries(fhout[0], kc), "writing unique boundaries file");
        _ACTION(save_kc(fhout[2], kc), "writing unique keycounts file");
    }
err:
    return res;
}

int
fa_index(struct seqb2_t* seq)
{
    struct gzfh_t* fhio[3] = { seq->fh, seq->fh + 1, seq->fh + 3};
    int len, res = -ENOMEM;
    char file[768];
    kct_t kc = {0};
    kc.ext = seq->readlength - KEY_WIDTH;
    unsigned mode;

    C char* ext[7] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub", ".uq"};

    if (fhio[0]->name) {
        len = strlen(fhio[0]->name);
    } else {
        ASSERT(seq->fh[2].name != NULL, return -EFAULT);
        fhio[0]->name = file;

        len = strlen(seq->fh[2].name);
        strncpy(fhio[0]->name, seq->fh[2].name, ++len);
        _ACTION0(reopen(fhio[0], ext[0], ext[5]), "")
    }
    ASSERT(len < 256, return -EFAULT, "filename too long: %s", fhio[0]->name);
    for (int i=1; i != 3; ++i) {
        if (fhio[i]->name == NULL)
            fhio[i]->name = &file[len*i];

        strncpy(fhio[i]->name, fhio[0]->name, len);
    }

    kc.ndxkct = _buf_init_arr_err(kc.ndxkct, KEYNT_BUFSZ_SHFT, return -ENOMEM);
    // first check whether unique boundary is ok.
    if (fhio[0]->fp) {
        mode = 2;
    } else {
         _ACTION0(reopen(fhio[0], ext[5], ext[4]), "trying to open %s instead", ext[4])
        mode = fhio[0]->fp != NULL;
    }
    if (mode) {
        _ACTION(load_boundaries(fhio[0], &kc), "loading boundary file %s", fhio[0]->name)
    }

    // keycount file did not exist, check twobit file.
    for (int i=0; i != 3; ++i) {
        _ACTION0(reopen(fhio[i], ext[i || mode == 2 ? 5 : 4], ext[i+1]), 
                    "%s does%s exist", ext[i+1], fhio[i]->fp ? "" : " not")
    }
    if (mode && fhio[0]->fp && fhio[1]->fp && fhio[2]->fp) {
         _ACTION(load_seqb2(fhio[0], &kc), "loading twobit sequence file")
         _ACTION(load_nextnts(fhio[1], &kc), "loading next Nts file")
         _ACTION(load_kc(fhio[2], &kc), "loading keycounts file")
    } else {
        bool found = false;
        for (int i=0; i != 3; ++i) {
            if (fhio[i]->fp) {
                found = true;
                EPR("%s file present, refusing to overwrite.", fhio[i]->name);
                rclose(fhio[i]);
            }
        }
        if (mode) {
            EPR("%s file present, refusing to overwrite.", ext[3+mode]);
            found = true;
        }
        if (found) goto err;
        EPR("starting from scratch.");
        _ACTION(fa_read(seq, &kc), "reading fasta")
    }
    if (mode < 2) {
        _ACTION(reopen(fhio[0], ext[1], ext[5]), "")
        _ACTION(reopen(fhio[2], ext[3], ext[6]), "")
        _ACTION(extd_uniqbnd(&kc, fhio), "extending unique boundaries")
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


