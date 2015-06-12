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

void
show_list(kct_t* kc, std::list<uint32_t> &bnd)
{
    std::list<uint32_t>::iterator it = bnd.begin();
    unsigned i = 0;
    for (it = bnd.begin(); it != bnd.end(); ++it) {
        Bnd* bd = &kc->bd[*it];
        EPR0("[%u (%u)]:\t%u\t+%u\t", i++, *it, bd->s, bd->l);
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

static int
decr_excise(kct_t const *const kc, const unsigned left)
{
    Walker* wlkr = kc->wlkr;
    uint64_t* wbuf = kc->wbuf;
    unsigned i = kc->ext;
    while (i-- > left) {

        // current bit pos and pending offsets are stored in walker.
        ASSERT(wbuf[i] != ~0ul, return -EFAULT, "%u/%u", i, left);
        dbg = wbuf[i] == dbgndx ?  dbg | 8 : dbg & ~8;
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
        EPQ0(dbg > 5, "before excision:");print_dna(*q, dbg > 5);
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
            EPQ0(dbg > 6, "became; next:");print_2dna(*q,q[1], dbg > 6);
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
    wbuf[i] = ~0ul;
    return left;
}

static unsigned iter = 0;

// if we cannot extend a range, update temporarily stored offsets.
static inline int update_wlkr(kct_t* kc, int left)
{
    uint64_t ndx = kc->wbuf[left];
    ASSERT(ndx < kc->kct_l, return -EFAULT, "[%u/%u]:%lx", left, kc->ext, ndx)
    ASSERT(kc->wlkr[ndx].tmp_count > 0, return -EFAULT);
    ++kc->wlkr[ndx].count;
    --kc->wlkr[ndx].tmp_count;
    kc->wbuf[left] = ~0ul;
    return 0;
}

static inline void merge(Bnd *dest, Bnd *next)
{
    EPQ0(dbg > 4, "Extending %u(-%u)", dest->s, dest->s + dest->l);
    dest->l = next->s - dest->s + next->l;
    dest->corr += next->corr;
    dest->dna = next->dna;
    EPQ(dbg > 4," to %u", dest->s + dest->l);
}

// XXX: inferiority ok? XXX: store deviant bit when unique,
// reverse seq could be required for mapping - maybe 
static int
ext_uq_bnd(kct_t* kc, Hdr* h, uint32_t lastx)
{
    // N-stretches, later also skips, future: SNVs, splice sites?
    int uqct = 0;
    // for realpos corrections (N-stretch or new contig) infior is zero.
    // otherwise the boundary describes a unique region - we can skip some.

    const char* hdr = kc->id + h->part[0];

    // one extra must be available for inter.
    _buf_grow0(kc->bd, 2ul);
    Bnd *inter = kc->bd + kc->bd_l; // running storage for boundary [ins or mod].
    Bnd *next = kc->bd + *kc->bdit;

    Bnd * last = kc->bd + lastx;
    uint32_t b2pos = last->s + last->l;

    uint64_t dna = last->dna; // first seq after skip
    uint64_t wx, ndx, rc = revcmp(dna);
    _get_ndx(ndx, wx, dna, rc);
    uint64_t infior = kc->kctndx[ndx] & R_INFIOR;
    if (infior == 0) infior += INFERIORITY;

    *inter = {.at_dna = dna, .dna = 0, .s = b2pos, .l = 0, .corr = last->corr};
    uint8_t b2;

    // boundary is considered as a first unique,
    // if we find a 2nd, we really want to extend this region
    const int ext = kc->ext;
    int left = ext;
    kc->wbuf[left-1] = ndx;
    ++kc->wlkr[ndx].tmp_count;
    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"",hdr,
            b2pos, next->s, next->s + next->l);
    print_dna(dna, dbg > 3);
    show_list(kc, h->bnd);
    //dbg = 7 * !strncmp(hdr, "GL000207.1", strlen(hdr));

    // for genomic pos add last->l.
    for (; b2pos < next->s;++b2pos, dna = _seq_next(b2, dna, rc)) { // until next event
        // this always stores the central Nt state regardless of uniqueness
        // wx gets strand and inferiority of key.
        ndx = _update_kctndx(kc, ndx, wx, dna, rc);
        unsigned ct = kc->kct[ndx]; // ts.idx of this key minus ts.idx of former is this keys' ts-length.
        uint64_t offs = ndx ? kc->kct[ndx-1] : 0;
        Walker* w = kc->wlkr + ndx;
        ct -= offs + w->excise_ct;

        offs += w->count + w->tmp_count;

        EPQ(dbg > 6, "offs:%lu\tend:%lu\tnext Nts (%lu part, %ux, byte %lu(%x)ndx:0x%lx):",
                offs, kc->kct[ndx], 4 - (offs & 3), ct, offs >> 2, kc->ts[offs >> 2], ndx);
        print_2dna(kc->ts[offs >> 2], kc->ts[offs >> 2] >> ((offs & 3) << 1), dbg > 6);
        EPQ0(dbg > 5, "%u:\t", b2pos); print_2dna(dna, rc, dbg > 5);

        // read next Nt from string
        b2 = (kc->ts[offs >> 2] >> ((offs & 3) << 1)) & 3;

        _verify_seq(b2pos, h, kc, "[%u]: sb2:%c, got %c, kctndx 0x%lx, [%u, %u]",
                return print_2dna(dna, rc), b2pos,b6(sb2<<1), b6(b2<<1), ndx, w->count, w->tmp_count)
        if (left) { // within uniq range: keep filling buffer
            if (--left) {
                if (left == ext - 1) {
                    EPQ(dbg > 6, "setting uniq jump destination");
                    // Sucessively overwritten until region completed - left became 0.
                    inter->dna = dna;
                    inter->l = b2pos - inter->s;
                }
                ASSERT(kc->wbuf[left-1] == ~0ul, return -EFAULT, "[%u/%u]", left-1,ext);
                kc->wbuf[left-1] = ndx;
                if (ct > 1ul) {
                    if ((infior + INFERIORITY) > wx) {
                         // A non-unique keys' inferiority must be gt neighbouring
                         // unique keys. Strand bit included, o well..
                         wx ^= infior + INFERIORITY;
                         kc->kctndx[ndx] ^= wx;
                     }
                    ++w->tmp_count;
                } else {
                    EPQ(dbg > 4, "2nd+ unique, extend range at %u, infior %u", b2pos, infior >> INFIOR_SHFT);
                    if (infior <= wx)
                         infior = wx & R_INFIOR; // only set strand bit below
                    // A uniqs' inferior must be ge neighbouring uniqs' infior.
                    wx ^= infior;
                    kc->kctndx[ndx] ^= wx;
                    w->count = 0;
                    w->tmp_count = 0;
                    ++uqct;
                    left = decr_excise(kc, left);
                    ASSERT(left >= 0, return -EFAULT, "left? %d", left);
                    if ((last->s + ext) >= inter->s) {
                        // join - may be undesirable for certain boundary types in the future
                        ASSERT(last != inter, return -EFAULT);
                        merge(last, inter);
                        *inter = {.at_dna = 0, .dna = 0, .s = last->s,
                            .l = 0, .corr = last->corr};
                    }
                    left = ext;
	            kc->wbuf[left-1] = ndx;
                    ++kc->wlkr[ndx].tmp_count;
                }
                continue;
            } else { //buf just became full, without extensions.
                infior = INFERIORITY;
                for (left = 0; left != ext; ++left) {
                    if (update_wlkr(kc, left) < 0)
                        return -EFAULT;
                }
                left = 0;
            }
        }
        if (ct > 1ul) {
            ++w->count;
        } else {
            EPQ(dbg > 3, "[%u]\t%u - %u\t", uqct, inter->s, inter->s + inter->l);

            EPQ(dbg > 5, "1st unique at %u", b2pos);
            w->count = 0;
            w->tmp_count = 0;
            inter->at_dna = dna;
            left = ext;
	    kc->wbuf[left-1] = ndx;
            ++kc->wlkr[ndx].tmp_count;
        }
    }
    ASSERT(b2pos == next->s, return -EFAULT, "%u, %u", b2pos, next->s);
    if (left--) { // now also last uniq must be included.
        inter->l = next->s - inter->s;
        ASSERT((inter->s + inter->l + ext) >= next->s, return -EFAULT, "%u %u", inter->s + inter->l, next->s);
        EPQ(dbg > 3, "[%u]\t%u - %u...\t", uqct, inter->s, next->s + next->l);
        print_2dna(last->dna, inter->at_dna, dbg > 1);

        left = decr_excise(kc, left);
        if (left < 0)
            return left;
        // if last and next have now become adjoining, merge
        EPQ(dbg > 3, "joining %u(-%u) til %u\t", inter->s,
            inter->s + inter->l, next->s + next->l);
        ASSERT(inter != next, return -EFAULT);
        merge(inter, next);
        kc->bdit = h->bnd.erase(kc->bdit);
        ASSERT(*--kc->bdit < kc->bd_l, return -EFAULT, "%lu >= %lu?", *kc->bdit, kc->bd_l);
        next = &kc->bd[*kc->bdit];
    }
    EPR(":::END:::");
    show_list(kc, h->bnd);
    return uqct;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    int uqct = 0;
    EPR("Processing %s (%lu)", kc->id + h->part[0],
            kc->bd[h->bnd.back()].s + kc->bd[h->bnd.back()].l);
    kc->bdit = h->bnd.begin();
    uint32_t corr = kc->bd[*kc->bdit].corr;

    for (uint32_t lastx = *kc->bdit; ++kc->bdit != h->bnd.end(); lastx = *kc->bdit) {
        int ret = ext_uq_bnd(kc, h, lastx);
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
    //dbg = 7;
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


