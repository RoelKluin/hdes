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
show_list(kct_t* kc, std::list<uint32_t> &bnd)
{
    std::list<uint32_t>::iterator it = bnd.begin();
    unsigned i = 0;
    for (it = bnd.begin(); it != bnd.end(); ++it) {
        Bnd* bd = &kc->bd[*it];
        EPR0("[%u (%u)]:\t%u\t+%u\t(+%u)\t", i++, *it, bd->s, bd->l, bd->corr);
        print_2dna(bd->at_dna, bd->dna);
    }
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
    _buf_free(kc->bd);
    _buf_free(kc->id);
    _buf_free(kc->ts);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->ndxkct);
}

static inline unsigned
pot_region_start(Bnd **C reg, C uint64_t dna, C uint32_t b2pos, C unsigned ext)
{
    reg[0]->s = b2pos;
    reg[0]->at_dna = dna;
    // A first uniq marks a potential start of a region
    reg[0]->dna = dna;
    reg[0]->l = 0;
    return b2pos + ext - max(reg[1]->s + reg[1]->l + ext, b2pos - ext);
}
static inline void
elongate_region(Bnd *C inter, C uint64_t dna, C uint32_t b2pos)
{
    inter->l = b2pos - inter->s;
    inter->dna = dna;
}
static inline int
end_region(kct_t *C kc, Bnd **C reg, Hdr* h)
{
    EPQ(dbg > 3 && reg[0]->l, "[%u]\t%u - %u\t", kc->bd_l, reg[0]->s, reg[0]->s + reg[0]->l);
    h->mapable += reg[0]->l;
    if (reg[0] != reg[1]) {
        reg[0]->corr = reg[1]->corr;
        _buf_grow0(kc->bd, 2ul);
        h->bnd.insert(kc->bdit, kc->bd_l);
        reg[1] = kc->bd + kc->bd_l++;
        reg[2] = kc->bd + *kc->bdit;
    }
    reg[0] = kc->bd + kc->bd_l; // no longer last at least
    reg[0]->l = reg[0]->s = 0;
    return 0;
}
static inline int
adjoining_boundary(kct_t *C kc, Bnd **C reg, Hdr *C h, C uint64_t dna, C uint32_t b2pos)
{
    reg[0]->dna = dna;
    reg[0]->l = b2pos - reg[0]->s;
    ASSERT(reg[0]->s + reg[0]->l == reg[2]->s, return -EFAULT);
    if (*reg == reg[1]) {
        EPQ(dbg > 4, "Removing boundary after loop [%u]", b2pos);
        // TODO: return *kc->bdit to a pool of to be inserted boundaries
        reg[1]->l += reg[2]->l;
        reg[1]->dna = reg[2]->dna;
        reg[1]->corr = reg[2]->corr;
        kc->bdit = h->bnd.erase(kc->bdit);
        reg[2] = kc->bd + *--kc->bdit;
    } else {
        reg[0]->corr = reg[1]->corr;
        EPQ (dbg > 4, "Next boundary joined after loop [%u]", b2pos);
        reg[2]->at_dna = reg[0]->at_dna;
        reg[2]->l += reg[2]->s - reg[0]->s;
        reg[2]->s = reg[0]->s;
        reg[2]->corr += reg[0]->corr;
    }
    return 0;
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
get_nextnt(kct_t C*C kc, C uint64_t nt)
{
    ASSERT(nt < (kc->ts_l << 2), return -EFAULT);
    return _GET_NEXT_NT(kc, nt); // next Nt
}

static inline int
mark_leaving(kct_t C*C kc, uint64_t *C k)
{
    int res = -EFAULT;
    if (PENDING(k)) {// not excised nor already handled
        *k -= ONE_PENDING - 1ul;
        uint64_t t = NEXT_NT_NR(k) + TSO_NT(k); // position of nextnt, to compare with former nextnt.
        ASSERT(t < (kc->ts_l << 2), goto err);

        EPQ(TSO_NT(k) == dbgtsoffs, "\t%c", b6(_GET_NEXT_NT(kc, t)<<1));

        if (IS_FIRST(k)) {
            if (IS_DISTINCT(k)) {
                EPQ(TSO_NT(k) == dbgtsoffs, " ==> distinct cleared");
                k[1] &= ~DISTINCT; // clear from last iteration
            } else if (kc->iter) {
                //if ok, then we can remove the kct iteration post loop?
                ASSERT(k[1] & MARKED, goto err, "MARKED not obsolete");
            }
        } else {
            ASSERT(t, goto err, "nt offset -1 (below)?");
            if (_GET_NEXT_NT(kc, t) != _GET_NEXT_NT(kc, t - 1)) {
                EPQ(TSO_NT(k) == dbgtsoffs, " ==> dbgtsoffs got distinct");
                k[1] |= DISTINCT;
            }
        }
    } else {
        EPQ(TSO_NT(k) == dbgtsoffs, "none pending or already handled");
        EPQ(TSO_NT(k) == dbgtsoffs, "*%s%s%s%s%sREMAIN:%lu\tNEXT_NT_NR:%lu\tPENDING:%lu\t",
            IS_FIRST(k) ? "FIRST\t" : "", ALL_SAME_NTS(k) ? "ALL_SAME\t" : "",
            IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",
            IS_LAST(k) ? "LAST\t" : "",REMAIN(k), NEXT_NT_NR(k), PENDING(k));
    }
    _MARK_LAST(k)
    res = 0;
err:
    if (res || TSO_NT(k) == dbgtsoffs || dbgtsoffs == -2ul) {
        EPR("<%s%s%s%sREMAIN:%lu(-1?)\tNEXT_NT_NR:%lu\tPENDING:%lu",
                ALL_SAME_NTS(k) ? "ALL_SAME\t" : "",
                IS_DISTINCT(k) ? "DISTINCT\t" : "", IS_UQ(k) ? "UQ\t" : "",
                IS_LAST(k) ? "LAST\t" : "",REMAIN(k), NEXT_NT_NR(k), PENDING(k));
        EPQ(res, "[%lu]", TSO_NT(k));
    }
    return res;
}
#define mark_leaving(kc, k) ({\
    EPQ0(TSO_NT(k) == dbgtsoffs, "%s +%u mark_leaving()", __FILE__, __LINE__);\
    mark_leaving(kc, k);\
})

static int
mark_all(kct_t C*C kc, unsigned i, unsigned const end)
{
    int res = 0;
    while (i != end)
        _EVAL(mark_leaving(kc, kc->kct_scope[i++]));
err:
    EPQ(res, "error during mark_all()");
    return res;
}


static int
decr_excise(kct_t C*C kc, uint64_t *C k)
{
    uint64_t nt, offs;
    uint8_t *q, *qe, c, t;
    int err = -EFAULT;

    ASSERT(IS_UQ(k) == false, goto out);

    if (ALL_SAME_NTS(k) == false) // only when previously incremented
        *k -= ONE_PENDING;

    k[1] -= ONE_CT; // one less remains.
    if (IS_LAST(k)) {//  src of movement would be target.
        err = 0;
        goto out;
    }

    offs = REMAIN(k), nt = TSO_NT(k);

    offs += nt; // remaining + ts offset => target of nt movement

    nt += NEXT_NT_NR(k); // add current position => src Nt (which is moved)

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
                IS_LAST(k) ? "LAST\t" : "",REMAIN(k), NEXT_NT_NR(k), PENDING(k));
        EPQ(err, "[%lu]", TSO_NT(k));
    }
    return err;

}
#define decr_excise(kc, k) ({\
    EPQ(TSO_NT(k) == dbgtsoffs, "%s +%u decr_excise()",\
            __FILE__, __LINE__);\
    decr_excise(kc, k);\
})

static int
decr_excise_all(kct_t C*C kc, uint64_t * k, uint64_t C*C kct, unsigned const end)
{
    int res = 0;
    uint64_t infior = max(*k, *kct);
    for (unsigned i = 1; i != end; ++i) {
        k = kc->kct_scope[i];
        if (k != kct) // .. don't elevate key if it is the one that became uniq
            update_max_infior(k, infior);
                _EVAL(decr_excise(kc, k));
        _MARK_LAST(k)
    }
err:
    EPQ(res, "before uniq %lu", TSO_NT(kct));
    return res;
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
ext_uq_bnd(kct_t *C kc, Hdr *C h, C uint32_t lastx)
{

    C char* hdr = kc->id + h->part[0];

    _buf_grow0(kc->bd, 2ul);    // one extra must be available for inter.
    Bnd *reg[3] = { kc->bd + lastx, kc->bd + lastx, kc->bd + *kc->bdit };

    uint64_t p, dna = reg[1]->dna;   // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t dummy[2] = {0ul, 0ul};
    uint64_t* k, *kct = dummy;
    unsigned i, index, ext = kc->ext - 1;
    uint32_t b2pos = reg[1]->s + reg[1]->l;
    int res;
    for (index = 0; index != ext; ++index)
        kc->kct_scope[index] = dummy;
    index = 0;

    //dbg = strncmp(hdr, "GL000207.1", strlen(hdr)) ? 3 : 5;
    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"", 
            hdr, b2pos, reg[2]->s, reg[2]->s + reg[2]->l);
    print_dna(dna, dbg > 3);

    while(b2pos < reg[2]->s) { // until next uniq region,  stretch or contig

        uint64_t nt;
        p = b2pos + h->s_s;
        /* process new key */
        _get_new_kct(kc, kct, dna, rc);
        // move in next key if last key indicates pending nextnt extension.
        /*k = kc->kct_scope[KC_LEFT_ROT_NEXT(ext, rot)];
        if (k != dummy && ALL_SAME_NTS(k) && IS_DISTINCT(k) == false) {
            EPR0("[b2pos:%u, tsoffs:%lu\t%lu], \t", b2pos, k[1] & B2POS_MASK, kct[1] & B2POS_MASK);
            print_2dna(dna, rc);
            EPR0("X");
            _EVAL(decr_excise(kc, kct, NULL, 0));
            ++*k;
            if (WAS_LAST(k)) { // all are now moved
                EPR0("Z");
                k[1] |= DISTINCT;
                mark_same_kct(kc, kct);
            }
            EPR0("\n");
            //KC_ROT(ext, rot); // not here.
        } */

        if (IS_UQ(kct)/* || NEXT_NT_NR(kct) == REMAIN(kct)*/) {

            *kct ^= (*kct ^ (p + 1)) & B2POS_MASK;
            _EVAL(get_nextnt(kc, kct[1]));
            nt = res;

            k = kc->kct_scope[0];
            if (IS_UQ(k) == false) { // then: first uniq
                //EPQ(TSO_NT(kct) == dbgtsoffs, "(kct == dbg: first uniq)");

                h->mapable += pot_region_start(reg, dna, b2pos, ext);
	        _EVAL(mark_all(kc, 0, index));
            } else {
                //EPQ(TSO_NT(kct) == dbgtsoffs, "(kct == dbg: 2nd or later uniq within scope)");
                // 2nd or later uniq within scope
                elongate_region(*reg, dna, b2pos);
		if (k != dummy)
        	    ++kc->uqct;

                // elevate all inbetween non-uniques by the max of these:
                _EVAL(decr_excise_all(kc, k, kct, index));
            }
	    index = 0;
        } else {

            if (ALL_SAME_NTS(kct)) {

                nt = _GET_NEXT_NT(kc, kct[1]);
                EPQ(TSO_NT(kct) == dbgtsoffs, "got %ul", nt);

            } else {
                _EVAL(get_nextnt(kc, NT_OFFS(kct)));
                nt = res;
                // NEXT_NT_NR must be instantly updated or next_nt won't match
                // twobit Nts for the same keys within scope.
                *kct += ONE_PENDING;
            }

	    if (index == ext) {
		i = 0;
		if (IS_UQ(kc->kct_scope[i])) {// first out of scope
		    _EVAL(end_region(kc, reg, h));
		    ++i;
		}
		_EVAL(mark_all(kc, i, index));
		index = 0;
	    }
        }
        _EVAL(check_nt(kc, nt, p & B2POS_MASK));

        ASSERT(nt != -1ul, res = -EFAULT; goto err, "ndxkct 0x%lx [%u]", *kct, b2pos);
        ++b2pos;
        kc->kct_scope[index++] = kct;

        dna = _seq_next(nt, dna, rc);
    }
    ASSERT(b2pos == reg[2]->s, res = -EFAULT; goto err, "[%u] %u", b2pos, reg[2]->s);
    k = kc->kct_scope[0];
    if (IS_UQ(k)) { // there was a unique in scope
        EPQ (dbg > 4, "Post loop boundary handling [%u]", b2pos);
        _EVAL(decr_excise_all(kc, k, kct, index));
	if (k != dummy)
            ++kc->uqct;
        _EVAL(adjoining_boundary(kc, reg, h, dna, b2pos));
    } else {
	_EVAL(mark_all(kc, 0, index));
        h->mapable -= max(reg[1]->s + reg[1]->l + ext, b2pos) - b2pos;
    }

    res = 0;
err:
    if (res < -1 || dbgtsoffs == -2ul) {
        EPR0("[b2pos:%u, tsoffs:%lu], \t", b2pos, kct[1] & B2POS_MASK);
        print_2dna(dna, rc);
    }
    return res;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    uint32_t lastx;
    EPR("Processing %s", kc->id + h->part[0]);
    kc->bdit = h->bnd.begin();
    h->mapable = 0u;

    for (lastx = *kc->bdit; ++kc->bdit != h->bnd.end(); lastx = *kc->bdit) {
        int ret = ext_uq_bnd(kc, h, lastx);
        if (ret < 0) return ret;
    }
    h->mapable += kc->bd[lastx].l; // last mapable region
    lastx = kc->bd[lastx].s + kc->bd[lastx].l;
#ifdef DEBUG
    if (dbg > 4)
        show_list(kc, h->bnd);
#endif
    EPQ(dbg > 2, "%s: %u/%u => %.2f%% mapable",
            kc->id + h->part[0], h->mapable, lastx,
            h->end_pos ? 100.0f * h->mapable / lastx : nanf("NAN"));

    return lastx;
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
    EPQ(dbg > 0, "extended %u unique ranges in iteration %u, using %u boundaries\n"
            "\t%lu/%lu => %.2f%% mapable. (%u pending)", kc->uqct, ++kc->iter, kc->bd_l,
            mapable, totNts, totNts ? 100.0f * mapable / totNts : nanf("NAN"), kc->pending);
    //dbg = 7;
    // FIXME: reset upon last occurance for non-uniq in inner loop
    for (uint64_t *k = kc->kct; k != kc->kct + kc->kct_l; k+=2) {
        if (IS_UQ(k) == false) {// at least one left (if unique position genomic 2bit)
#if defined(TEST_CLEARANCE)
            ASSERT(NEXT_NT_NR(k) == 0ul, return -EFAULT, "[%lu]", TSO_NT(k));
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


