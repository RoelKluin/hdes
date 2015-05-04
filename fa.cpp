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
#include <string.h> // memset()
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
        free(h->s);
        free(h->part);
        delete h;
    }
    _buf_free(kc->bd);
    _buf_free(kc->id);
    _buf_free(kc->ts);
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
        EPQ0(dbg > 2, "before excision:");print_dna(*q, dbg > 2);
        uint8_t c = ((nt & 3) + 1) << 1;
        uint8_t t = *q & ((1ul << c) - 1ul); // capture this and past Nts.
        if (c) c -= 2;
        *q = ((*q ^ t) >> 2) | (*q & ((1ul << c) - 1ul)); // excise out Nt.
        t >>= c;
        EPQ(dbg > 2, "moving %c from %lu to %lu", b6(t << 1), nt, offs);

        // mask to cover this and past nucleotide.
        EPQ(dbg > 2, "=====> excised <======, became %x",
                *q | (q != qe ? (q[1] & 3) << 6 : 0));
        wbuf[i] = ~0ul;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            EPQ0(dbg > 2, "became; next:");print_2dna(*q,q[1], dbg > 2);
            *++q >>= 2;
        }
        // append excised Nt to end.
        offs = (offs & 3) << 1;        // set shift
        EPQ(dbg > 2, "%lu, %c", offs, b6(t << 1));
        t <<= offs;                    // move excised in position
        offs = *q & ((1u << offs) - 1u); // below 2bits were shifted correctly
        *q = ((*q ^ offs) << 2) ^ t ^ offs;  // move top part back up, add rest.
        EPQ0(dbg > 2, "after append:");print_dna(*q, dbg > 2);
    }
    wbuf[i] = ~0ul;
    return left;
}

static unsigned iter = 0;

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
    const char* dbgtid = "GL000210.1 ";//"GL000207.1";//MT"; //GL000239.1";
    const char* hdr = kc->id + h->part[0];
    //dbg = (strncmp(hdr, dbgtid, strlen(dbgtid)) == 0) * 2;

    uint64_t dna = last->dna; // first seq after skip
    uint64_t rc = revcmp(dna);
    uint32_t this_infior;
    uint64_t ndx, wx;
    _update_kctndx(kc, ndx, wx, this_infior, dna, rc);

    uint32_t infior = max(last->i + 1, this_infior); //XXX

    // one extra must be available for inter.
    _buf_grow0(kc->bd, 2ul);
    Bnd *inter = &kc->bd[kc->bd_l];
    *inter = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = infior};

    // boundary is considered as a first unique,
    // if we find a 2nd, we really want to extend this region
    const int ext = kc->ext;
    int left = ext;
    EPQ0(dbg > 1, "----[\t%s%s:(%u-%u)..(%u-%u)\t@%u\t]----\tdna:",
            strlen(hdr) < 8 ? "\t":"",hdr, last->s, last->s+last->l, next->s,
            next->s+next->l, pos);
    print_dna(dna, dbg > 1);

    for (;pos < next->s;++pos) { // until next event
        EPQ(dbg & 4, "observed dbgndx at %u", pos);
        EPQ0(dbg > 2, "%u:\t", pos);
        print_2dna(dna, rc, dbg > 2);

        Walker* w = kc->wlkr + wx;
        if (left) { // within uniq range: keep filling buffer
            if (left == ext) {
                //EPR("position after 2nd unique, or each extended.");
                // what we'll jump to, in subsequent iterations.
                // This is inserted (once) when region is complete (left becomes 0) -
                // otherwise it's overwritten successively.
                inter->dna = dna;
                inter->l = pos - inter->s;
            }
            if (--left) {
                ASSERT(kc->wbuf[left] == ~0ul, return -EFAULT, "[%u]", pos & print_dna(dna));
                kc->wbuf[left] = wx;
                if (infior > this_infior) {
                    // replace this inferiority
                    kc->kctndx[ndx] &= INDEX_MASK;
                    kc->kctndx[ndx] |= (uint64_t)infior << (STRAND_SHFT + 1);
                }
            } else { // actually buf just became full, without extension(s not already handled)..
                EPQ(inter->s > 0xfffffff, "%u", -inter->s);
                infior = 0;
                // add tmp_count for each - we cannot extend it in this iteration
                for (left = ext; --left;) {
                    ++kc->wlkr[kc->wbuf[left]].count;
                    --kc->wlkr[kc->wbuf[left]].tmp_count;
                    kc->wbuf[left] = ~0ul;
                }

                if (inter->l != 0) {
                    // we did insert a range
                    EPQ0(dbg > 1, "\n[%u-%u]\t", inter->s, inter->s + inter->l - 1);print_2dna(inter->at_dna, inter->dna, dbg > 1);
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
        }
        uint64_t offs = wx != 0 ? kc->kct[wx-1] : 0;
        unsigned ct = kc->kct[wx] - offs - w->excise_ct;
        offs += w->count + w->tmp_count;

        EPQ(dbg > 2, "offs:%lu\tend:%lu\tnext Nts (%lu part, %ux, byte %lu(%x)ndx:0x%lx):",
                offs, kc->kct[wx], 4 - (offs & 3), ct, offs >> 2, kc->ts[offs >> 2], ndx);
        print_2dna(kc->ts[offs >> 2], kc->ts[offs >> 2] >> ((offs & 3) << 1), dbg > 2);

        // read next Nt from string
        uint8_t b2 = (kc->ts[offs >> 2] >> ((offs & 3) << 1)) & 3;

        if (h->s) { // compare to sequence string - not stored on disk.
            uint8_t sb2 = (h->s[pos>>2] >> ((pos & 3) << 1)) & 3;
            ASSERT(b2 == sb2, return print_2dna(dna, rc),
                    "[%u]: expected %c, got %c for ndx 0x%lx, [%u, %u]",
                    pos,b6(sb2<<1), b6(b2<<1), ndx, w->count, w->tmp_count);
        }
        dna = _seq_next(b2, dna, rc);
        // this always stores the central Nt state regardless of uniqueness
        _update_kctndx(kc, ndx, wx, this_infior, dna, rc);
        if (ct > 1ul) {
            if (left == 0) ++w->count;
            else ++w->tmp_count;
            continue;
        }
        EPQ0(dbg > 1, "%u,", pos);
        w->count = 0;
        w->tmp_count = 0;

        ++infior;
        if (left == 0) { EPQ(dbg > 2, "this is a first unique.");
            *inter = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = max(infior, this_infior + 1u)};
            left = ext;
            infior = inter->i + 1;
            continue;
        }
        // insert a new range.
        ++uqct;
        left = decr_excise(kc, ext, left);
        if (left < 0) return left;
        left = ext;
    }
    ASSERT(pos == next->s, return -EFAULT);
    EPQ(dbg > 2, "---");
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

int
fn_convert(struct gzfh_t* fhio, const char* search, const char* replace)
{
    char* f = strstr(fhio->name, search);
    if (f == NULL) return 0;
    strncpy(f, replace, strlen(replace) + 1);
    return 1;
}

int fa_index(struct seqb2_t* seq, uint32_t blocksize)
{
    struct gzfh_t* fhio = seq->fh;
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    size_t t;
    int res = -ENOMEM, ret = -1;
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];
    kct_t kc;
    kc.ext = seq->readlength - KEY_WIDTH;
    kc.kctndx = _buf_init_arr(kc.kctndx, KEYNT_BUFSZ_SHFT);

    const char* ndxact[4] = {".fa", "_x1.gz", "_x2.gz"};

    if (fhio->name != NULL && ((res = strlen(fhio->name)) > 255)) {
        strncpy(file, fhio->name, res);
    } else {
        res = strlen(fhin->name);
        ASSERT(strncmp(fhin->name, "stdin", res) != 0, return -1);
        strcpy(file, fhin->name);
    }
    fhio->name = &file[0];

    // first check whether last file exists
    if (!fn_convert(fhio, ndxact[0], ndxact[2]))
        return -1;

    fhio->fp = fopen(fhio->name, "r");
    if (fhio->fp == NULL) {

        // last file did not exist, check 2nd last.
        if (!fn_convert(fhio, ndxact[2], ndxact[1]))
            return -1;

        fhio->fp = fopen(fhio->name, "r");

        if (fhio->fp == NULL) { // write both first output file and last.
            res = set_io_fh(fhio, blocksize, 0);
            ASSERT(res >= 0, goto out);

            kc.kct = _buf_init_err(kc.kct, 16, goto out);
            kc.kce = _buf_init_err(kc.kce, 8, goto out);
            kc.bd = _buf_init_err(kc.bd, 1, goto out);

            kc.id = _buf_init_err(kc.id, 5, goto out);

            memset(kc.kctndx, ~0ul, KEYNT_BUFSZ * sizeof(*kc.kctndx));

            if (is_gzfile) {
                g = fhin->io;
                gc = (int (*)(void*))&gzgetc;
                ungc = (int (*)(int, void*))&gzungetc;
            } else {
                g = fhin->fp;
                gc = (int (*)(void*))&fgetc;
                ungc = (int (*)(int, void*))&ungetc;
            }
            /* TODO: load dbSNP and known sites, and mark them. */
            res = fa_kc(&kc, g, gc, ungc);
            EPR("done reading and intializing keycounts");
            ASSERT(res >= 0, goto out);

            res = kct_convert(&kc);
            EPR("done converting keycounts");
            _buf_free(kc.kce);
            ASSERT(res >= 0, goto out);

            res = write1(fhio, &kc);
            EPR("done writing(1) %s to disk", fhio->name);
        } else { // read first file, write last.
            res = set_io_fh(fhio, blocksize, 2);
            ASSERT(res >= 0, goto out);
            res = restore1(fhio, &kc);
            EPR("done loading(1) %s from disk", fhio->name);
        }
        ASSERT(res >= 0, goto out);
        res = rclose(fhio);
        ASSERT(res >= 0, goto out);

        t = kc.kct_l * sizeof(Walker);
        kc.wlkr = (Walker*)malloc(t);
        ASSERT(kc.wlkr != NULL, goto out)
        memset(kc.wlkr, 0ul, t);

        t = kc.ext * sizeof(uint64_t);
        kc.wbuf = (uint64_t*)malloc(t);
        if (kc.wbuf == NULL) {
            free(kc.wlkr);
            goto out;
        }
        memset(kc.wbuf, ~0ul, t);

        do { // until no no more new uniques
            res = ext_uq_iter(&kc);
        } while (res > 0);
        _buf_free(kc.wlkr);
        _buf_free(kc.wbuf);
        ASSERT(res >= 0, goto out);

        // reopen for writing processed boundaries (2nd file)
        if (!fn_convert(fhio, ndxact[1], ndxact[2]))
            return -1;
        res = set_io_fh(fhio, blocksize, 0);
        ASSERT(res >= 0, goto out, "failed to open %s", fhio->name);

        res = write1(fhio, &kc);
        EPR("done writing(2) %s to disk", fhio->name);

    } else {
        // read last file
        res = set_io_fh(fhio, blocksize, 2);
        ASSERT(res >= 0, goto out);

        res = restore1(fhio, &kc);
        EPR("done loading(2) %s from disk", fhio->name);
    }
    ret = res != 0;
out:
    EPQ(ret, "an error occured:%d", res);
    free_kc(&kc);
    return ret;

}


