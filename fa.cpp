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

static unsigned iter = 0;
static uint64_t maxinfior = 0;

static inline void
update_max_infior(uint64_t* kct, uint64_t infior)
{
    uint64_t k = *kct;
    if (infior > k) {
        maxinfior = max(maxinfior, infior);
        *kct ^= (k ^ infior) & INFIOR_MASK;
    }
}

static inline uint64_t
get_nextnt(kct_t const *const kc, uint64_t* kct, uint64_t p)
{
    uint64_t t = IS_UQ(kct) ? 0 : (*kct)++;         // passed this key, if not yet unique

    *kct ^= -!t & (*kct ^ p) & STRAND_POS; // otherwise set b2pos
    t = (t + kct[1]) & B2POS_MASK;         // next Nt offset for this key

    t = (kc->ts[t>>2] >> ((t&3) << 1)) & 3; // next Nt
    if (kc->s == NULL)
        return t;

    p = --p & B2POS_MASK; // zero based.
    p = (kc->s[p>>2] >> ((p&3) << 1)) & 3;
    if (t == p)
        return t;

    WARN("nextnt didn't match twobit at b2pos %lu: sb2:%c, got %c, ndxkct 0x%lx",
            p, b6(p<<1), b6(t<<1), kct - kc->kct);
    return -1ul;
}


static int
decr_excise(kct_t const *const kc, uint64_t* kct, uint64_t infior)
{
    uint64_t nt = kct[1] -= ONE_CT; // excise it: one less remains..
    uint64_t offs = nt >> BIG_SHFT;
    nt &= B2POS_MASK; // ts offset
    offs += nt;
    //update_max_infior(kct, infior);

    nt += --(*kct) & B2POS_MASK; // ..instead of passing
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


// reverse seq could be required for mapping - maybe 
// N-stretches, later also skips, future: SNVs, splice sites?
static int
ext_uq_bnd(kct_t* kc, Hdr* h, uint32_t lastx)
{

    const char* hdr = kc->id + h->part[0];

    _buf_grow0(kc->bd, 2ul);    // one extra must be available for inter.
    Bnd *next = kc->bd + *kc->bdit;
    Bnd *last = kc->bd + lastx;
    Bnd *inter = last;          // running storage for boundary [ins or mod].

    uint32_t b2pos = last->s + last->l;
    uint64_t dna = last->dna;   // first seq after skip
    uint64_t rc = revcmp(dna);
    running r = {0};

    unsigned rest = kc->ext;
    uint64_t ndx, t = _ndxkct_and_infior(ndx, t, dna, rc);
    uint64_t* kct;

    //else dbg = strncmp(hdr, "GL000207.1", strlen(hdr)) ? 3 : 5;
    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"", 
            hdr, b2pos, next->s, next->s + next->l);
    print_dna(dna, dbg > 3);

    while(b2pos < next->s) { // until next uniq region,  stretch or contig

        // insertion of kct index for new complement insensitive index
        kct = kc->kct + kc->ndxkct[ndx];

        t |= ++b2pos + h->s_s; // position in sam is one based.
        if (IS_UQ(kct)) {
            r.infior = t;
        } else { // for non-unique ensure infior is above last uniq infior.
            update_max_infior(kct, r.infior + INFIOR);
        }

        // TODO: if all nextNts are the same increase keylength.

        t = get_nextnt(kc, kct, t);
        if (t == -1ul)
            return -EFAULT & print_2dna(dna, rc);

        // within scope of unique key or when leaving it
        switch(rest) {
case 1:     EPQ(dbg > 3 && inter->l, "[%u]\t%u - %u\t", kc->bd_l, inter->s, inter->s + inter->l);
            h->mapable += inter->l;
            if (inter != last) {
                inter->corr = last->corr;
                _buf_grow0(kc->bd, 2ul);
                h->bnd.insert(kc->bdit, kc->bd_l);
                last = kc->bd + kc->bd_l++;
                next = kc->bd + *kc->bdit;
            }
            inter = kc->bd + kc->bd_l; // no longer last at least
            inter->l = inter->s = 0;
            --rest;
            break;
case 0:     r.last = r.rot;
            break;
default:    --rest;
            r.rot |= -!r.rot & kc->ext; // if zero, rotate to kc->ext
            kc->wbuf[--r.rot] = kct;
            if (IS_UQ(kct)) {
                unsigned j = r.last;
                j |= -!j & kc->ext; // skip first (unique).
                while(--j != r.rot) {
                    uint64_t* k = kc->wbuf[j];
                    // no movement if only one left, write genomic position
                    if (IS_UQ(k) == false) {
                        decr_excise(kc, k, r.infior + INFIOR);
                    } else {
                        uint64_t tp = b2pos + h->s_s + j - r.rot;
                        if (j > r.rot) tp -= kc->ext;
                        if ((*k & B2POS_MASK) != tp) {
                            EPQ((*k & B2POS_MASK) != 1ul,"%lu, %lu", *k & B2POS_MASK, tp);
                            *k ^= (*k ^ tp) & B2POS_MASK;
                        }
                    }
                    j |= -!j & kc->ext;
                }
            }
        }
        // when in scope of an uniq key - keylength minus keywidth distance,
        if (IS_UQ(kct)) {
            if (rest) {
                // 2nd or later uniq within scope => region insertion or update.
                if (inter != last)
                    ++kc->uqct;
                // Sucessively overwritten until region completed - left became 0.
                inter->dna = dna;
                inter->l = b2pos - inter->s - 1;
            } else {
                //XXX: make sure this is not off by one or two.
                h->mapable += b2pos - max(last->s + last->l + kc->ext, b2pos - kc->ext);
                h->mapable += kc->ext;
                inter->s = b2pos - 1;
                inter->at_dna = dna;
                // A first uniq marks a potential start of a region
                inter->dna = dna;
                inter->l = 0;
            }
            r.last = r.rot % kc->ext;
            rest = kc->ext;
        } else if (r.rot == r.last) {
            r.infior = 0;
        }
        dna = _seq_next(t, dna, rc);
        t = _ndxkct_and_infior(ndx, t, dna, rc);

    }
    ASSERT(b2pos == next->s, return -EFAULT, "%u, %u", b2pos, next->s);
    if (rest) {
        EPQ (dbg > 4, "Post loop boundary handling at %u", b2pos);
        if (r.rot != r.last) {
            ++kc->uqct;

            for(unsigned j = r.last; j != r.rot;) {
                j |= -!j & kc->ext;
                kct = kc->wbuf[--j];
                // no movement if only one left, write genomic position
                if (IS_UQ(kct) == false) {
                    decr_excise(kc, kct, r.infior + INFIOR);
                } else {
                    uint64_t tp = b2pos + h->s_s + j - r.rot;
                    if (j > r.rot) tp -= kc->ext;
                    if ((*kct & B2POS_MASK) != tp) {
                        EPQ((*kct & B2POS_MASK) != 1ul,"%lu, %lu", *kct & B2POS_MASK, tp);
                        *kct ^= (*kct ^ tp) & B2POS_MASK;
                    }
                }
            }
        }
        inter->dna = dna;
        inter->l = b2pos - inter->s;
        ASSERT(inter->s + inter->l == next->s, return -EFAULT)
        if (inter == last) {
            EPQ(dbg > 4, "Removing boundary after loop at %u", b2pos);
            // TODO: return *kc->bdit to a pool of to be inserted boundaries
            last->l += next->l;
            last->dna = next->dna;
            last->corr = next->corr;
            kc->bdit = h->bnd.erase(kc->bdit);
            next = kc->bd + *--kc->bdit;
        } else {
            inter->corr = last->corr;
            EPQ (dbg > 4, "Next boundary joined after loop at %u", b2pos);
            next->at_dna = inter->at_dna;
            next->l += next->s - inter->s;
            next->s = inter->s;
            next->corr += inter->corr;
        }
    } else {
        h->mapable -= max(last->s + last->l + kc->ext, b2pos) - b2pos;
    }
    return 0;
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
    kc->wbuf = (uint64_t**)malloc(t * sizeof(uint64_t*));
    ASSERT(kc->wbuf != NULL, goto err);

    do { // until no no more new uniques
        res = ext_uq_iter(kc);
    } while (res > 0);
    _buf_free(kc->wbuf);
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


