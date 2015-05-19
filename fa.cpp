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
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"

static void
show_list(kct_t* kc, std::list<uint32_t> &bnd)
{
    std::list<uint32_t>::iterator it = bnd.begin();
    unsigned i = 0;
    for (it = bnd.begin(); it != bnd.end(); ++it) {
        Bnd* bd = &kc->bd[*it];
        EPR0("[%u (%u)]:\t%u-%u\t%u\t", i++, *it, bd->s, bd->s + bd->l, bd->l);
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

static void
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

static int
decr_excise(kct_t const *const kc, unsigned i, const unsigned left)
{
    Walker* wlkr = kc->wlkr;
    uint64_t* wbuf = kc->wbuf;
    while (--i > left) {

        // current bit pos and pending offsets are stored in walker.
        ASSERT(wbuf[i] != ~0ul, return -EFAULT, "%u/%u", i, left);
        dbg = wbuf[i] == dbgndx ?  dbg | 4 : dbg & ~4;
        Walker* w = wlkr + wbuf[i];
        ASSERT(w->tmp_count > 0, return -EFAULT);

        // get ts offset
        uint64_t* kct = kc->kct;
        uint64_t nt = wbuf[i] != 0 ? kct[wbuf[i] - 1] : 0;
        nt += w->count + --w->tmp_count;
        uint64_t offs = kct[wbuf[i]] - ++w->excise_ct;

        // excise out twobits, i.e. shift 2bits above current 2bit down.
        uint8_t* q = &kc->ts[nt >> 2];
        uint8_t* qe = &kc->ts[offs >> 2];

        // shift for past Nt;
        EPQ0(dbg > 3, "before excision:");print_dna(*q, dbg > 3);
        uint8_t c = ((nt & 3) + 1) << 1;
        uint8_t t = *q & ((1ul << c) - 1ul); // capture this and past Nts.
        if (c) c -= 2;
        *q = ((*q ^ t) >> 2) | (*q & ((1ul << c) - 1ul)); // excise out Nt.
        t >>= c;
        EPQ(dbg > 3, "moving %c from %lu to %lu", b6(t << 1), nt, offs);

        // mask to cover this and past nucleotide.
        EPQ(dbg > 3, "=====> excised <======, became %x",
                *q | (q != qe ? (q[1] & 3) << 6 : 0));
        wbuf[i] = ~0ul;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            EPQ0(dbg > 3, "became; next:");print_2dna(*q,q[1], dbg > 3);
            *++q >>= 2;
        }
        // append excised Nt to end.
        offs = (offs & 3) << 1;        // set shift
        EPQ(dbg > 2, "%lu, %c", offs, b6(t << 1));
        t <<= offs;                    // move excised in position
        offs = *q & ((1u << offs) - 1u); // below 2bits were shifted correctly
        *q = ((*q ^ offs) << 2) ^ t ^ offs;  // move top part back up, add rest.
        EPQ0(dbg > 3, "after append:");print_dna(*q, dbg > 3);
    }
    wbuf[i] = ~0ul;
    return left;
}

static unsigned iter = 0;

// if we cannot extend a range, update temporarily stored offsets.
static inline int update_wlkr(kct_t* kc, int left)
{
    uint64_t ndx = kc->wbuf[left];
    ASSERT(kc->wlkr[ndx].tmp_count > 0, return -EFAULT);
    ++kc->wlkr[ndx].count;
    --kc->wlkr[ndx].tmp_count;
    kc->wbuf[left] = ~0ul;
    return 0;
}

static inline void merge(Bnd *dest, Bnd *next)
{
    EPQ(dbg > 1, "Extending %u-%u til %u", dest->s, next->s, next->s + next->l);
    dest->l += next->l;
    dest->dna = next->dna;
}

// XXX: inferiority ok? XXX: store deviant bit when unique,
// reverse seq could be required for mapping - maybe 
static int
ext_uq_bnd(kct_t* kc, Hdr* h, Bnd *last)
{
    // N-stretches, later also skips, future: SNVs, splice sites?
    int uqct = 0;
    uint32_t pos = last->s + last->l;
    Bnd *next = &kc->bd[*kc->bdit];
    const char* hdr = kc->id + h->part[0];

    uint64_t dna = last->dna; // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t infior = last->i & ~INDEX_MASK;
    // one extra must be available for inter.
    _buf_grow0(kc->bd, 2ul);
    Bnd *inter = &kc->bd[kc->bd_l]; // running storage.
    *inter = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = infior};
    uint8_t b2;

    // boundary is considered as a first unique,
    // if we find a 2nd, we really want to extend this region
    const int ext = kc->ext;
    int left = ext;
    EPQ0(dbg > 3, "----[\t%s%s:(%u-%u)..(%u-%u)\t@%u\t]----\tdna:",
            strlen(hdr) < 8 ? "\t":"",hdr, last->s, last->s+last->l, next->s,
            next->s+next->l, pos);
    print_dna(dna, dbg > 3);

    for (;pos < next->s;++pos, dna = _seq_next(b2, dna, rc)) { // until next event
        // this always stores the central Nt state regardless of uniqueness
        uint64_t ndx, wx; // wx gets strand and inferiority of key.
        ndx = _update_kctndx(kc, ndx, wx, dna, rc);
        unsigned ct = kc->kct[ndx];
        uint64_t offs = ndx ? kc->kct[ndx-1] : 0;
        Walker* w = kc->wlkr + ndx;
        ct -= offs + w->excise_ct;

        offs += w->count + w->tmp_count;

        EPQ(dbg > 3, "offs:%lu\tend:%lu\tnext Nts (%lu part, %ux, byte %lu(%x)ndx:0x%lx):",
                offs, kc->kct[ndx], 4 - (offs & 3), ct, offs >> 2, kc->ts[offs >> 2], ndx);
        print_2dna(kc->ts[offs >> 2], kc->ts[offs >> 2] >> ((offs & 3) << 1), dbg > 3);

        // read next Nt from string
        EPQ(dbg & 4, "observed dbgndx at %u", pos);
        EPQ0(dbg > 3, "%u:\t", pos);
        print_2dna(dna, rc, dbg > 3);

        b2 = (kc->ts[offs >> 2] >> ((offs & 3) << 1)) & 3;

        if (kc->s) { // compare to sequence string - not stored on disk.
            uint64_t p = h->s_s + pos;
            uint8_t sb2 = (kc->s[p>>2] >> ((p & 3) << 1)) & 3;
            ASSERT(b2 == sb2, return print_2dna(dna, rc),
                    "[%u]: expected %c, got %c for ndx 0x%lx, [%u, %u]",
                    pos,b6(sb2<<1), b6(b2<<1), ndx, w->count, w->tmp_count);
        }
        if (left) { // within uniq range: keep filling buffer
            if (--left) {
                if (left == ext - 1) {
                    //EPR("position after 2nd unique, or each extended.");
                    // what we'll jump to, in subsequent iterations.
                    // This is inserted (once) when region is complete (left becomes 0) -
                    // otherwise it's overwritten successively.
                    inter->dna = dna;
                    inter->l = pos - inter->s;
                }
                ASSERT(kc->wbuf[left] == ~0ul, return -EFAULT, "[%u]", pos & print_dna(dna));
                kc->wbuf[left] = ndx;
                if (ct > 1ul) {
                    if ((infior + INFERIORITY) > wx) {
                         // A non-unique keys' inferiority must be gt neighbouring
                         // unique keys. Strand bit included, o well..
                         wx ^= infior + INFERIORITY;
                         kc->kctndx[ndx] ^= wx;
                     }
                    ++w->tmp_count;
                } else { // another unique, extend range
                    if (infior > wx) {
                         // A unique keys' inferiority must be ge neighbouring
                         wx ^= infior; // unique keys, so update.
                     } else {          // otherwise update inferiority,
                         infior = wx & ~(1ul << (STRAND_SHFT + 1));
                         wx ^= infior; // and only set strand bit below
                     }
                    kc->kctndx[ndx] ^= wx;
                    EPQ0(dbg > 2, "%u,", pos);
                    w->count = 0;
                    w->tmp_count = 0;
                    ++uqct;
                    left = decr_excise(kc, ext, left);
                    if (left < 0) return left;
                    left = ext;
                }
                continue;
            }
            // else buf just became full, without extensions.
            EPQ(inter->s > 0xfffffff, "%u", -inter->s);
            infior = 0;
            for (left = ext; --left;) {
                if (update_wlkr(kc, left) < 0)
                    return -EFAULT;
            }

            if (inter->l != 0) {
                // we did insert a range
                EPQ0(dbg > 2, "\n[%u-%u]\t", inter->s, inter->s + inter->l - 1);print_2dna(inter->at_dna, inter->dna, dbg > 1);
                //show_list(kc, h->bnd);
                if ((last->s + last->l) < inter->s) {
                    h->bnd.insert(kc->bdit, kc->bd_l);
                    _buf_grow0(kc->bd, 2ul);
                    last = kc->bd + kc->bd_l++;
                    inter = last + 1;
                } else { // join - may be undesirable for certain boundary types in the future
                    merge(last, inter);
                }
                next = &kc->bd[*kc->bdit];
                inter->l = 0;
            }
        }
        if (ct > 1ul) {
            ++w->count;
        } else {
            infior = wx & ~(1ul << (STRAND_SHFT + 1));
            wx ^= infior; // and only set strand bit below

            EPQ0(dbg > 2, "%u,", pos);
            w->count = 0;
            w->tmp_count = 0;
            *inter = {.at_dna = 0, .dna = 0, .s = pos, .l = 0, .i = infior};
            left = ext;
        }
    }
    ASSERT(pos == next->s, return -EFAULT);
    if (left) {
        ASSERT((inter->s + inter->l + ext) >= next->s, return -EFAULT);
        EPQ0(dbg > 1, "\nExtending %u-%u til %u\t(next:%u-%u)\t", last->s, last->s+last->l, inter->s+inter->l, next->s, next->s+next->l); print_2dna(last->dna, inter->at_dna, dbg > 1);
        next->l += next->s - inter->s;
        next->s = inter->s;
        next->at_dna = inter->at_dna;
        //XXX: why decrement left?
        left = decr_excise(kc, ext, --left);
        if (left < 0)
            return left;
        // if last and next have now become adjoining, merge the two
        // XXX why not (last->s + last->l + ext > next->s) ? (assertion error later, but why?)
        if ((last->s + last->l >= next->s)) {
            merge(last, next);
            kc->bdit = h->bnd.erase(kc->bdit);
            --kc->bdit;
        }
    }
    return uqct;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    int uqct = 0;
    EPR("Processing %s", kc->id + h->part[0]);
    kc->bdit = h->bnd.begin();

    for (Bnd *last = &kc->bd[*kc->bdit]; ++kc->bdit != h->bnd.end(); last = &kc->bd[*kc->bdit]) {
        int ret = ext_uq_bnd(kc, h, last);
        if (ret < 0) return ret;
        uqct += ret;
    }
#ifdef DEBUG
    for (unsigned t = 0; t != kc->kct_l; ++t) {
        Walker* w = kc->wlkr + t;
        ASSERT(w->tmp_count == 0u, return -EFAULT, "%u/%lu", t, kc->kct_l);
    }
    if (dbg > 1)
        show_list(kc, h->bnd);
#endif

    return uqct;
}

static int
ext_uq_iter(kct_t* kc)
{
    int uqct = 0;
    for (std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
        uqct += ret;
    }
    EPR("extended %u unique ranges in iteration %u", uqct, iter++);
    for (Walker *w = kc->wlkr; w != &kc->wlkr[kc->kct_l]; ++w) {
        ASSERT(w->tmp_count == 0u, return -EFAULT);
        w->count = 0u;
    }
    return uqct;
}


static int
extd_uniqbnd(kct_t* kc, struct gzfh_t* fhout)
{
    int res = -ENOMEM;
    size_t t = kc->kct_l;
    kc->wlkr = (Walker*)malloc(t * sizeof(Walker));
    ASSERT(kc->wlkr != NULL, return res);
    while (t--) kc->wlkr[t] = {0u, 0u, 0u};

    t = kc->ext;
    kc->wbuf = (uint64_t*)malloc(t * sizeof(uint64_t));
    ASSERT(kc->wbuf != NULL, goto err);

    while (t--) kc->wbuf[t] = ~0ul;

    do { // until no no more new uniques
        res = ext_uq_iter(kc);
    } while (res > 0);
    _buf_free(kc->wbuf);
    _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
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

    kc.kctndx = _buf_init_arr(kc.kctndx, KEYNT_BUFSZ_SHFT);
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


