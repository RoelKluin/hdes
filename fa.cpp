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
#include <errno.h> // ENOMEM
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
    for (Kct* k = kc->kct; k != &kc->kct[kc->kct_l]; ++k)
        if (k->seq.m > 3)
            free(k->p.b2);

    std::map<char*, Hdr*, Tid>::iterator it;
    for (it = kc->hdr.begin(); it != kc->hdr.end(); ++it) {
        free(it->second->s);
        free(it->second->part);
        delete it->second;
    }
    _buf_free(kc->bd);
    _buf_free(kc->id);
    _buf_free(kc->kct);
    _buf_free(kc->kce);
    _buf_free(kc->wlkr);
    _buf_free(kc->wbuf);
    _buf_free(kc->kcsndx);
}

static const unsigned long dbgndx = UNINITIALIZED;

static int
decr_excise(kct_t const *const kc, unsigned i, const unsigned left)
{
    Walker* wlkr = kc->wlkr;
    uint32_t* wbuf = kc->wbuf;
    while (--i > left) {

        // current bit pos and pending offsets are stored in walker.
        ASSERT(wbuf[i] != UNINITIALIZED, return -EFAULT, "%u/%u", i, left);
        Kct *x = &_get_kct(kc, wbuf[i]);
        Walker* w = &get_w(wlkr, kc, wbuf[i]);
        // excise out twobits, i.e. shift 2bits above current 2bit
        // down. current bit pos is stored in walker.
        uint8_t* q, *qe;
        ASSERT(w->tmp_count > 0, return -EFAULT);
        --w->tmp_count;
        uint64_t t = w->count;
        if (x->seq.m == 3) {
            q = &x->seq.b2[t >> 2];
            qe = &x->seq.b2[(x->seq.l-1) >> 2];
            --x->seq.l;
        } else {
            q = &x->p.b2[t >> 2];
            qe = &x->p.b2[(x->p.l-1) >> 2];
            // could convert to Kct.seq, but why bother?
            --x->p.l;
        }
        // mask to cover this and past nucleotide.
        t = (1 << (((w->count & 3) + 1) << 1)) - 1;
        *q = ((*q & (t ^ 0xff)) >> 2) | (*q & (t >> 2));
        EPQ(wbuf[i] == dbgndx, "=====> excised <======, became %x",
                *q | (q != qe ? (q[1] & 3) << 6 : 0));
        wbuf[i] = UNINITIALIZED;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            *++q >>= 2;
        }
    }
    //Walker* w = &get_w(wlkr, kc, wbuf[i]);
    //ASSERT(w->tmp_count > 0, return -EFAULT);
    //--w->tmp_count;
    //EPQ(wbuf[i] == dbgndx, "=====> after excised <======");
    wbuf[i] = UNINITIALIZED;

    while (i != 0) {
        if (wbuf[i] != UNINITIALIZED) {
            Walker* w = &get_w(wlkr, kc, wbuf[i]);
            ASSERT(wbuf[i] == UNINITIALIZED, return print_ndx(wbuf[i]), "index %u/%u:%u", i, left, w->tmp_count);
        }
        --i;
    }
    return left;
}

static int dbg = 1;
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

    uint64_t t, dna = last->dna; // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t ndx = _getxtdndx(kc, ndx, dna, rc);
    Walker* w = &get_w(kc->wlkr, kc, ndx);
    uint32_t infior = max(last->i + 1, w->infior); //XXX
    if (dbg > 2)
        show_list(kc, h->bnd);

    // one extra must be available for inter.
    _buf_grow0(kc->bd, 2ul);
    Bnd *inter = &kc->bd[kc->bd_l];
    *inter = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = infior};

    // boundary is considered as a first unique,
    // if we find a 2nd, we really want to extend this region
    const int ext = kc->ext;
    int left = ext;
    EPQ0(dbg > 0, "----[\t%s%s:(%u-%u)..(%u-%u)\t@%u\t]----\tdna:",
            strlen(hdr) < 8 ? "\t":"",hdr, last->s, last->s+last->l, next->s,
            next->s+next->l, pos);
    print_dna(dna, dbg > 0);

    for (;pos < next->s;++pos) { // until next event
        EPQ(ndx == dbgndx, "observed dbgndx at %u", pos);
        ASSERT(_getxtdndx0(kc, ndx) < kc->kct_l, return -print_dna(dna), "at %u, 0x%lx", pos, ndx);
        EPQ0(dbg > 2, "%u:\t", pos);
        print_dna(dna, dbg > 2);

        w = &get_w(kc->wlkr, kc, ndx);
        if (left) { // keep filling the buffer
            if (left == ext) {
                //EPR("position after 2nd unique, or each extended.");
                // what we'll jump to, in subsequent iterations.
                // This is inserted (once) when region is complete (left becomes 0) -
                // otherwise it's overwritten successively.
                inter->dna = dna;
                inter->l = pos - inter->s;
            }
            if (--left) {
                ASSERT(kc->wbuf[left] == UNINITIALIZED, return -EFAULT, "[%u]", pos & print_dna(dna));
                kc->wbuf[left] = ndx;
                if (infior > w->infior)
                    w->infior = infior;
            } else { // actually buf just became full, without extension(s not already handled)..
                EPQ(inter->s > 0xfffffff, "%u", -inter->s);
                infior = 0;
                // add tmp_count for each - we cannot extend it in this iteration
                for (left = ext; --left;) {
                    get_w(kc->wlkr, kc, kc->wbuf[left]).count++;
                    --get_w(kc->wlkr, kc, kc->wbuf[left]).tmp_count;
                    kc->wbuf[left] = UNINITIALIZED;
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
        t = w->count + w->tmp_count;
        Kct *y = &_get_kct(kc, ndx);
        uint8_t b2;
        if (y->seq.m == 3 || y->seq.m == 2) { //EPR("Kct in seq format");
            b2 = (y->seq.b2[t >> 2] >> ((t & 3) << 1)) & 3;
            t = y->seq.l;
        } else {
            b2 = (y->p.b2[t >> 2] >> ((t & 3) << 1)) & 3; // get next 2bit
            t = y->p.l; // can also be 1 (unique), unless we decide to convert back upon decrement
        }
        if (h->s) { // not set after loading from disk.
            uint8_t sb2 = (h->s[pos>>2] >> ((pos & 3) << 1)) & 3;
            ASSERT(b2 == sb2, return print_2dna(dna, rc), "[%u]: expected %c, got %c for ndx 0x%lx", pos,b6(sb2<<1), b6(b2<<1), ndx);
        }
        dna = _seq_next(b2, dna, rc);
        ndx = _getxtdndx(kc, ndx, dna, rc);
        if (t > 1ul) {
            if (left == 0) ++w->count;
            else ++w->tmp_count;
            continue;
        }
        EPQ0(dbg > 1, "%u,", pos);
        // found unique, store orientation of 2nd bit of central Nt.
        // if last bit is set, the key is in the wrong orientation.
        // in any case - before any revcmping, the central bit needs to be
        // put back in to generate the sequence at this site.
        y->seq.m ^= (y->seq.m & 1) ^ !(dna & KEYNT_STRAND);
        ++infior;
        if (left == 0) { EPQ(dbg > 2, "this is a first unique.");
            *inter = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = max(infior, w->infior + 1u)};
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
    //ASSERT(dna == next->at_dna, return -print_2dna(dna, next->at_dna), "[%u] => enddna mismatch", pos);
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

    } else {
        EPQ(dbg > 1, "Not joining\n");
    }
//                pos += kc->bd[*(kc->bdit)].l; // add skipped Nts
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
    for (unsigned t = 0; t != kc->kct_l; ++t) { //XXX
        Walker* w = kc->wlkr + t;
        //XXX
        ASSERT(w->tmp_count == 0u, return -EFAULT, "%u/%u", t, kc->kct_l);
    }
    if (dbg > 1)
        show_list(kc, h->bnd);

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
    int res, ret = -1;
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];
    kct_t kc;
    kc.ext = seq->readlength - KEY_WIDTH;
    kc.kcsndx = _buf_init_arr(kc.kcsndx, KEYNT_BUFSZ_SHFT);

    const char* ndxact[4] = {".fa", "_x1.gz"};

    if (fhio->name != NULL && ((res = strlen(fhio->name)) > 255)) {
        strncpy(file, fhio->name, res);
    } else {
        res = strlen(fhin->name);
        ASSERT(strncmp(fhin->name, "stdin", res) != 0, return -1);
        strcpy(file, fhin->name);
    }
    fhio->name = &file[0];
    if (!fn_convert(fhio, ndxact[0], ndxact[1]))
        return -1; // skip if this file was not provided on the commandline

    fprintf(stderr, "== %s(%lu)\n", fhio->name, fhio - seq->fh);
    res = set_io_fh(fhio, blocksize, 0);
    if (res < 0) { fprintf(stderr, "== set_io_fh failed:%d", res); goto out; }

    if (res == 1) {
        kc.kct = _buf_init(kc.kct, 16);
        kc.kce = _buf_init(kc.kce, 8);
        kc.bd = _buf_init(kc.bd, 1);

        // FIXME: assigning 1 Mb for reference id, to keep realloc from happening.
        // (realloc would invalidate pointers inserted in kc.hdr)
        kc.id = _buf_init(kc.id, 20);

        memset(kc.kcsndx, UNINITIALIZED, KEYNT_BUFSZ * sizeof(kc.kcsndx[0]));

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
        EPR("left fa_kc");fflush(NULL);
        ASSERT(res >= 0, goto out);

        _buf_free(kc.kce);

        // write to disk, TODO: load from disk.
        res = write1(fhio, &kc);
        //ASSERT(res >= 0, goto out);
        //res = fhio->close(fhio->io);
    } else {
        res = restore1(fhio, &kc);
        ASSERT(res >= 0, goto out);
        res = fhio->close(fhio->io);
    }
    ASSERT(res >= 0, goto out);

    t = kc.kct_l * sizeof(Walker);
    kc.wlkr = (Walker*)malloc(t);
    memset(kc.wlkr, 0ul, t);

    t = kc.ext * sizeof(uint32_t);
    kc.wbuf = (uint32_t*)malloc(t);
    memset(kc.wbuf, UNINITIALIZED, t);

    do { // until no no more new uniques
        res = ext_uq_iter(&kc);
    } while (res > 0);
    ASSERT(res >= 0, goto out);

    ret = res != 0;

out:
    free_kc(&kc);
    return ret;

}


