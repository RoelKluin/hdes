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

static inline int
print_ndx(uint64_t dna, bool dbg = true)
{
    if (dbg == false) return -1;
    uint64_t rc = dna & KEYNT_TRUNC_UPPER;
    dna ^= rc ^ (rc << 1) ^ KEYNT_STRAND;
    rc = revcmp(dna);
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('|', stderr);
    for (unsigned t = KEY_WIDTH; t--; rc >>= 2)
        fputc(b6((rc & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
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
    _buf_free(kc->kctndx);
}

//FIXME also keep infior above a uniq keys' infior after this key.
static int
decr_excise(kct_t const *const kc, const int left)
{
    Walker* wlkr = kc->wlkr;
    uint64_t* wbuf = kc->wbuf;
    unsigned keep_dbg = 0;
    for(int i = kc->ext; i --> left ;) {

        // current next_Nt pos and pending offsets are stored in walker.
        uint64_t ndx = wbuf[i] & INDEX_MASK;
        ASSERT(ndx != ~0ul, return -EFAULT, "%u/%u", i, kc->ext);
        dbg = ndx == dbgkctndx ?  dbg | 8 : dbg & ~8;
        keep_dbg |= dbg;
        Walker* Next_nts = wlkr + ndx;
        ASSERT(Next_nts->count > 0, return -EFAULT);

        // get ts offset
        uint64_t* kct = kc->kct;
        uint64_t nt = ndx != 0 ? kct[ndx - 1] : 0;
        nt += --Next_nts->count;
        uint64_t offs = kct[ndx] - ++Next_nts->excise_ct;

        // excise out twobits, i.e. shift 2bits above current 2bit down.
        uint8_t* q = &kc->ts[nt >> 2];
        uint8_t* qe = &kc->ts[offs >> 2];

        // shift for past Nt;
        EPQ0(dbg > 5, "before excision:");print_2dna(*q,q[1], dbg > 5);
        uint8_t c = ((nt & 3) + 1) << 1;
        uint8_t t = *q & ((1ul << c) - 1ul); // capture this and past Nts.
        if (c) c -= 2;
        *q = ((*q ^ t) >> 2) | (*q & ((1ul << c) - 1ul)); // excise out Nt.
        t >>= c;
        EPQ(dbg > 5, "moving %c from %lu to %lu", b6(t << 1), nt, offs);

        EPQ(dbg > 6, "=====> excised <======, became %x",
                *q | (q != qe ? (q[1] & 3) << 6 : 0));
        // mask to cover this and past nucleotide.
        wbuf[i] = ~0ul;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            EPQ0(dbg > 6, "became; next:");print_2dna(*q,q[1]>>2, dbg > 6);
            *++q >>= 2;
        }
        // append excised Nt to end.
        offs = (offs & 3) << 1;        // set shift
        EPQ(dbg > 5, "%lu, %c", offs, b6(t << 1));
        t <<= offs;                    // move excised in position
        offs = *q & ((1u << offs) - 1u); // below 2bits were shifted correctly
        *q = ((*q ^ offs) << 2) ^ t ^ offs;  // move top part back up, add rest.
        EPQ0(dbg > 5, "after append:");print_dna(*q, dbg > 5);
    }
    //if(left != kc->ext) wbuf[left] = ~0ul;
    return keep_dbg;
}

static unsigned iter = 0;

// if we cannot extend a range, update temporarily stored offsets.
static inline int update_wlkr(kct_t* kc, int left)
{
    uint64_t ndx = kc->wbuf[left] & INDEX_MASK;
    if (ndx == dbgkctndx) {
        EPR("-- cleared dbgkctndx 0x%lx --", dbgkctndx);
        dbg |= 8;
    }
    ASSERT(ndx < kc->kct_l, return -EFAULT, "[%u/%u]:%lx", left, kc->ext, ndx)
    kc->wbuf[left] = ~0ul;
    return 0;
}

/*
 * Update strand and next_nt offset, when in scope of an uniq key - keylength minus
 * keywidth distance, ct==2 - keep this keys' inferiority above that infior.
 * A first uniq marks a potential start of a region, a 2nd or later uniq, within
 * scope - keylength minus keywidth distance - triggers an region insertion or
 * update. Also retrieve for this index the offset stored in kctndx high bits.
 */
static inline int eval_ndx(kct_t* kc, running *r, uint64_t *kctndx,
        uint64_t wx, unsigned ct)
{
    EPQ(dbg > 6, "path ct:%u", ct);
    switch(ct) { //(ct == 1) | ((left > 0) << 1);
case 2: if ((r->infior + INFERIORITY) > wx) {
             // A non-unique keys' inferiority must be gt neighbouring
             // unique keys. Strand bit included, but that shouldn't matter.
             wx ^= r->infior + INFERIORITY;
             *kctndx ^= wx;
        }
        _store_ndx(kc, r->left, *kctndx);
        break;
case 3:{// A uniqs' inferior must be ge neighbouring uniqs' infior.
        int ret  = decr_excise(kc, r->left);
        ASSERT(ret >= 0, return -EFAULT, "left:%d", ret);
        dbg |= ret;
        EPQ (dbg > 7, "dbg excision");
        wx ^= r->infior;
        *kctndx ^= wx; // XXX: no -ge test before infior replacement?
        ++kc->uqct;

        if (r->infior <= wx)
            r->infior = wx & R_INFIOR; // only set strand bit below
       }
case 1: r->left = kc->ext;
        _store_ndx(kc, r->left, *kctndx);
    }
    return *kctndx;
}

// XXX: inferiority ok? XXX: store deviant bit when unique,
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
    unsigned ct = 1; // treat last boundary as a first unique
    running r = {0};

    if (dbg > 5) show_list(kc, h->bnd);
    //else dbg = strncmp(hdr, "GL000207.1", strlen(hdr)) ? 3 : 5;
    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"", 
            hdr, b2pos, next->s, next->s + next->l);
    print_dna(dna, dbg > 3);

    while(b2pos < next->s) { // until next event
        EPQ(dbg > 6, "next path ct:%u", ct);
        uint64_t ndx, t = _kctndx_and_infior(ndx, t, dna, rc);
        EPQ0(dbg > 5, "%u:\t", b2pos); print_2dna(dna, rc, dbg > 5);

        // get offset to first of keys' next Nt
        int kct_i = eval_ndx(kc, &r, kc->kctndx + ndx, t, ct);
        ASSERT(kct_i >= 0, return -EFAULT, "(%s:%u)", hdr, b2pos);
        t = kct_i ? kc->kct[kct_i-1] : 0;

        // ts.idx of this key minus ts.idx of former is this keys' next-Nts-length.
        Walker* Next_nts = kc->wlkr + kct_i;
        ct = kc->kct[kct_i] - t - Next_nts->excise_ct;
        ASSERT(ct != 0, return -EFAULT);
        ct = ct == 1ul; // if only one, the key has become uniq

        t += Next_nts->count++; // offset to current next Nt.

        EPQ(dbg > 6, "offs:%lu\tnext Nts (%lu part, %ux, byte %lu(%x)):",
            t, 4 - (t&3), ct, t >> 2, kc->ts[t>>2]);
        print_2dna(kc->ts[t>>2], kc->ts[t>>2] >> ((t&3) << 1), dbg > 6);

        t = (kc->ts[t>>2] >> ((t&3) << 1)) & 3;
        if (kc->s) {
            uint64_t p = b2pos + h->s_s;
            p = (kc->s[p>>2] >> ((p&3) << 1)) & 3;
            if (t != p) {
                WARN("assertion 'next_b2 != sb2' failed [%u]: sb2:%c, got %c, kctndx 0x%x",
                    b2pos,b6(p<<1), b6(t<<1), kct_i);
                return print_2dna(dna, rc);
            }
        }
        switch(r.left) {
case 0:     break; // otherwise within uniq range: keep filling buffer
case 1:     EPQ(dbg > 3 && inter->l, "[%lu]\t%u - %u\t", kc->bd_l, inter->s, inter->s + inter->l);
            if (inter != last) {
                h->mapable += inter->l;
                inter->corr = last->corr;
                _buf_grow0(kc->bd, 2ul);
                h->bnd.insert(kc->bdit, kc->bd_l++);
                last = kc->bd + lastx;
                next = kc->bd + *kc->bdit;
            }
            for (unsigned i = kc->ext; i--> 0 ;) {
                if (update_wlkr(kc, i) < 0)
                    return -EFAULT;
            }
            inter = kc->bd + kc->bd_l; // no longer last at least
            inter->l = inter->s = 0;
            r.infior = 0;
            --r.left;
            break;
default:    if (r.left-- == kc->ext) {
                EPQ(dbg > 6, "setting uniq jump destination");
                // Sucessively overwritten until region completed - left became 0.
                inter->dna = dna;
                inter->l = b2pos - inter->s;
            }
            ct |= 2;
        }
        dna = _seq_next(t, dna, rc);
        ++b2pos;
        if (ct == 1) {
            EPQ(dbg > 5, "1st unique at %u", b2pos);
            inter->s = b2pos;
            inter->at_dna = dna;
        }
    }
    h->mapable += last->l;
    EPQ(dbg > 5, "Last inter: s:%u, l:%u", inter->s, inter->l);
    ASSERT(b2pos == next->s, return -EFAULT, "%u, %u", b2pos, next->s);
    if (r.left || ct == 3) {
        EPQ (dbg > 4, "Post loop boundary handling at %u", b2pos);
        if (r.left) {
            ++kc->uqct;
            decr_excise(kc, r.left);
        }

        inter->dna = dna;
        inter->l = b2pos - inter->s;
        if (last->s + last->l == next->s) {
            EPQ(dbg > 4, "Removing boundary after loop at %u", b2pos);
            EPQ(dbg > 6, "last:%u +%u &%u, inter:%u +%u &%u next:%u +%u &%u", last->s, last->l, last->corr, inter->s, inter->l, inter->corr, next->s, next->l, next->corr);
            // TODO: return *kc->bdit to a pool of to be inserted boundaries
            h->mapable -= last->l; // it will again be added later.
            last->l += next->l;
            last->dna = next->dna;
            last->corr = next->corr;
            kc->bdit = h->bnd.erase(kc->bdit);
            next = kc->bd + *--kc->bdit;
            EPQ(dbg > 5, "next became:%u +%u &%u", next->s, next->l, next->corr);
        } else if (inter != last) {
            inter->corr = last->corr;
            if (inter->s + inter->l != next->s) {
                EPQ (dbg > 5, "Boundary insertion after loop at %u", b2pos);
                _buf_grow0(kc->bd, 2ul);
                h->bnd.insert(kc->bdit, kc->bd_l++);
            } else {
                EPQ (dbg > 4, "Next boundary joined after loop at %u", b2pos);
                next->at_dna = inter->at_dna;
                next->l += next->s - inter->s;
                next->s = inter->s;
                next->corr += inter->corr;
            }
        }
    }
    return 0;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    kc->uqct = 0u;
    uint32_t lastx = h->bnd.back();
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
    EPQ(dbg > 0, "extended %u unique ranges in iteration %u, using %lu boundaries\n"
            "\t%lu/%lu => %.2f%% mapable", kc->uqct, iter++, kc->bd_l, mapable, totNts,
            totNts ? 100.0f * mapable / totNts : nanf("NAN"));
    //dbg = 7;
    for (Walker *w = kc->wlkr; w != kc->wlkr + kc->kct_l; ++w)
        w->count = 0u;

    return kc->uqct;
}


static int
extd_uniqbnd(kct_t* kc, struct gzfh_t* fhout)
{
    int res = -ENOMEM;
    size_t t = kc->kct_l;
    kc->wlkr = (Walker*)malloc(t * sizeof(Walker));
    ASSERT(kc->wlkr != NULL, return res);
    while (t--) kc->wlkr[t] = {0u, 0u};

    t = kc->ext;
    kc->wbuf = (uint64_t*)malloc(t * sizeof(uint64_t));
    ASSERT(kc->wbuf != NULL, goto err);

    while (t--) kc->wbuf[t] = ~0ul;

    do { // until no no more new uniques
        res = ext_uq_iter(kc);
    } while (res > 0);
    _buf_free(kc->wbuf);
    if (res == 0) {
        _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
    }
err:
    _buf_free(kc->wlkr);
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

    const char* ext[6] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub"};

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

    kc.kctndx = _buf_init_arr_err(kc.kctndx, KEYNT_BUFSZ_SHFT, return -ENOMEM);
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
        _ACTION(extd_uniqbnd(&kc, fhio[0]), "extending unique boundaries")
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


