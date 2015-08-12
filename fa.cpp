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
pot_region_start(Bnd **reg, const uint64_t dna, const uint32_t b2pos, const unsigned ext)
{
    reg[0]->s = b2pos;
    reg[0]->at_dna = dna;
    // A first uniq marks a potential start of a region
    reg[0]->dna = dna;
    reg[0]->l = 0;
    return b2pos + ext - max(reg[1]->s + reg[1]->l + ext, b2pos - ext);
}
static inline void
elongate_region(Bnd *inter, uint64_t dna, uint32_t b2pos)
{
    inter->l = b2pos - inter->s;
    inter->dna = dna;
}
static inline int
end_region(kct_t *const kc, Bnd **reg, Hdr* h)
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
adjoining_boundary(kct_t *const kc, Bnd **reg, Hdr* h, const uint64_t dna, const uint32_t b2pos)
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

static unsigned iter = 0;
static uint64_t maxinfior = 0;

static inline void
update_max_infior(uint64_t* kct, uint64_t* exception, uint64_t infior)
{
    if (kct != exception) {// .. don't elevate key if it is the same as the one at r.rot.
        infior += INFIOR;
        if (infior > *kct) {
            maxinfior = max(maxinfior, infior);
            *kct ^= (*kct ^ infior) & INFIOR_MASK;
        }
    }
}

static inline uint64_t
get_nextnt(kct_t const *const kc, uint64_t nt, uint64_t p)
{
    ASSERT(nt < (kc->ts_l << 2), return -1ul);
    nt = (kc->ts[nt>>2] >> ((nt&3) << 1)) & 3; // next Nt
    if (kc->s == NULL)
        return nt;

    ASSERT(p < (kc->s_l << 2), return -1ul, "%lu/%lu", p, kc->s_l);
    p = (kc->s[p>>2] >> ((p&3) << 1)) & 3;
    if (nt == p)
        return nt;

    WARN("nextnt didn't match twobit: sb2:%c, got %c", b6(p<<1), b6(nt<<1));
    return -1ul;
}

static int
decr_excise(kct_t const *const kc, uint64_t* kct, uint64_t* exception, uint64_t infior)
{

    if (kct != exception) {// .. don't elevate key if it is the same as the one at r.rot.
        // update max inferiority
        infior += INFIOR;
        if (infior > *kct) {
            maxinfior = max(maxinfior, infior);
            *kct ^= (*kct ^ infior) & INFIOR_MASK;
        }
        // if exception and unique, the same key occured multiple times in same region
        if (IS_UQ(kct))
            return 0;
    }

    uint64_t nt = kct[1] -= ONE_CT;
    --*kct;
    uint64_t offs = nt >> BIG_SHFT; // remaining
    nt &= B2POS_MASK; // ts offset
    offs += nt; // add current position

    nt += *kct & B2POS_MASK;
    if (nt == offs) {
        // also if at last nextNt we can skip.
        // TODO: maybe check here whether all nextNts are the same - for
        // extended key in next iterations, set remain to 0? ts offset obsolete.
        return 0;
    }
    ASSERT((offs >> 2) < kc->ts_l, return -EFAULT);
    ASSERT((nt >> 2) + 1 < kc->ts_l, return -EFAULT);

    // excise out twobits, i.e. shift 2bits above current 2bit down.
    uint8_t* q = &kc->ts[nt >> 2];
    uint8_t* qe = &kc->ts[offs >> 2];
    ASSERT(q <= qe, return -EFAULT, "%lu\t%lu", nt, offs);

    // shift for past Nt;
    uint8_t c = ((nt & 3) + 1) << 1;
    uint8_t t = *q & ((1ul << c) - 1ul); // capture this and past Nts.
    if (c) c -= 2;
    *q = ((*q ^ t) >> 2) | (*q & ((1ul << c) - 1ul)); // excise out Nt.
    t >>= c;

    // mask to cover this and past nucleotide.
    while (q != qe) {
        *q |= (q[1] & 3) << 6; //XXX
        *++q >>= 2;
    }
    // append excised Nt to end.
    offs = (offs & 3) << 1;        // set shift
    t <<= offs;                    // move excised in position
    offs = *q & ((1u << offs) - 1u); // below 2bits were shifted correctly
    *q = ((*q ^ offs) << 2) ^ t ^ offs;  // move top part back up, add rest.
    return 0;
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
ext_uq_bnd(kct_t* kc, Hdr* h, uint32_t lastx)
{

    const char* hdr = kc->id + h->part[0];

    _buf_grow0(kc->bd, 2ul);    // one extra must be available for inter.
    Bnd *reg[3] = { kc->bd + lastx, kc->bd + lastx, kc->bd + *kc->bdit };

    uint64_t dna = reg[1]->dna;   // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t dummy[2] = {0ul, ONE_CT};
    running r = {0};
    uint64_t* kct;
    uint32_t b2pos = reg[1]->s + reg[1]->l;
    int res;
    kc->kct_scope[0] = dummy;

    //else dbg = strncmp(hdr, "GL000207.1", strlen(hdr)) ? 3 : 5;
    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"", 
            hdr, b2pos, reg[2]->s, reg[2]->s + reg[2]->l);
    print_dna(dna, dbg > 3);

    while(b2pos < reg[2]->s) { // until next uniq region,  stretch or contig

        /* process leaving key */
        KC_ROT(kc, r.rot);

        /* process new key */
        uint64_t p = _get_new_kct(kc, kct, b2pos + h->s_s, dna, rc, r.rot);

        uint64_t nt = kct[1];

        if (IS_UQ(kct)) {
            r.infior = *kct ^= (*kct ^ (p+1)) & STRAND_POS; //position in sam is one-based.
            // when the leaving and new key were both unique we cannot consider it
            // a complete region.
            if (r.last_uq == r.rot) {

                _ACTION(end_region(kc, reg, h), "");
                h->mapable += pot_region_start(reg, dna, b2pos, kc->ext);

            } else if (IS_UQ(kc->kct_scope[r.last_uq])) {
                // 2nd or later uniq within scope
                if (kc->kct_scope[r.last_uq] != dummy)
                    ++kc->uqct;
                elongate_region(*reg, dna, b2pos);

                // elevate all inbetween non-uniques by the max of these:
                uint64_t infior = max(*kc->kct_scope[r.last_uq], *kct);
                while (KC_ROT(kc, r.last_uq) != r.rot) {
                    _ACTION(decr_excise(kc, kc->kct_scope[r.last_uq], kct, infior), "");
                }
            } else {
                r.last_uq = r.rot;
                h->mapable += pot_region_start(reg, dna, b2pos, kc->ext);
            }
        } else {
            // ts_offs + passed => current Nt, passed this key
            nt += (*kct)++;
            if (IS_UQ(kc->kct_scope[r.last_uq])) {
                update_max_infior(kct, NULL, r.infior);
            } else if (r.last_uq == r.rot) {
                r.infior = 0;
                _ACTION(end_region(kc, reg, h), "");
            } else {
                r.last_uq = r.rot;
            }
        }
        nt = get_nextnt(kc, nt & B2POS_MASK, p & B2POS_MASK);
        ASSERT(nt != -1ul,  return print_2dna(dna, rc), "ndxkct 0x%lx [%u]", *kct, b2pos);
        ++b2pos;
        dna = _seq_next(nt, dna, rc);
    }
    ASSERT(b2pos == reg[2]->s, return -EFAULT, "[%u] %u", b2pos, reg[2]->s);

    if (IS_UQ(kc->kct_scope[r.last_uq])) { // there was a unique in scope
        EPQ (dbg > 4, "Post loop boundary handling [%u]", b2pos);

        if (r.rot != r.last_uq) {
            ++kc->uqct;
            uint64_t infior = kc->kct_scope[r.last_uq][1];

            while (r.last_uq != r.rot) {
                _ACTION(decr_excise(kc, kc->kct_scope[KC_ROT(kc, r.last_uq)], NULL, infior), "");
            }
        }
        _ACTION(adjoining_boundary(kc, reg, h, dna, b2pos), "");
    } else {
        h->mapable -= max(reg[1]->s + reg[1]->l + kc->ext, b2pos) - b2pos;
    }
    res = 0;
err:
    EPQ(res < -1, "[%u]", b2pos);
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
    kc->uqct = 0u;
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
            "\t%lu/%lu => %.2f%% mapable\n\tmaxinfior: %lx", kc->uqct, iter++, kc->bd_l,
            mapable, totNts, totNts ? 100.0f * mapable / totNts : nanf("NAN"), maxinfior);
    //dbg = 7;
    // FIXME: reset upon last occurance for non-uniq in inner loop
    for (uint64_t *k = kc->kct; k != kc->kct + kc->kct_l; k+=2)
        if (!IS_UQ(k)) // at least one left (if unique position genomic 2bit)
            *k &= REMAIN_MASK; // reset position tracking

    return kc->uqct;
}


static int
extd_uniqbnd(kct_t* kc, struct gzfh_t** fhout)
{
    int res = -ENOMEM;
    size_t t;

    // was ndx for storage, unset for pos, strand & infiority;
    for (unsigned i=0u; i != kc->kct_l; i += 2)
        kc->kct[i] = 0ul;

    t = kc->ext;
    kc->kct_scope = (uint64_t**)malloc(t * sizeof(uint64_t*));
    ASSERT(kc->kct_scope != NULL, goto err);

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

    const char* ext[7] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub", ".uq"};

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


