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

//FIXME also keep infior above a uniq keys' infior after this key.
static int
decr_excise(kct_t const *const kc, const unsigned left)
{
    Walker* wlkr = kc->wlkr;
    uint64_t* wbuf = kc->wbuf;
    unsigned i = kc->ext;
    while(i--> left) {

        // current next_Nt pos and pending offsets are stored in walker.
        uint64_t ndx = wbuf[i];
        ASSERT(ndx != ~0ul, return -EFAULT, "%u/%u", i, left);
        dbg = ndx == dbgkctndx ?  dbg | 8 : dbg & ~8;
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
    return left;
}

static unsigned iter = 0;

// if we cannot extend a range, update temporarily stored offsets.
static inline int update_wlkr(kct_t* kc, int left)
{
    uint64_t ndx = kc->wbuf[left];
    ASSERT(ndx < kc->kct_l, return -EFAULT, "[%u/%u]:%lx", left, kc->ext, ndx)
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

/*
 * Update strand and next_nt offset, if in scope of an uniq key - keylength minus
 * keywidth distance, ct==2 - keep this keys' inferiority above that infior.
 * A first uniq marks a potential start of a region, a 2nd or later uniq, within
 * scope - keylength minus keywidth distance - triggers an region insertion or
 * update. Also retrieve for this index the offset stored in kctndx high bits.
 */
static inline uint64_t eval_ndx(kct_t* kc, running *r, Bnd *inter,
        uint64_t dna, uint64_t rc, uint32_t b2pos)
{
    uint64_t wx, ndx;
    _get_ndx(ndx, wx, dna, rc);

    // store strand orientation and inferiority in wx
    wx <<= STRAND_SHFT - KEY_WIDTH;
    wx |= kc->kctndx[ndx] & ~INDEX_MASK;

    // get keycount index.
    ndx = kc->kctndx[ndx] & INDEX_MASK;

    Walker* Next_nts = kc->wlkr + ndx;
    unsigned add = Next_nts->count;

    EPQ(dbg > 6, "path ct:%u", r->ct);
    switch(r->ct) { //(ct == 1) | ((left > 0) << 1);
case 1: EPQ(dbg > 5, "1st unique at %u", b2pos);
        inter->s = b2pos;
        inter->at_dna = dna;
        r->left = kc->ext;
        _incr_at_ndx(kc, r->left, ndx);
        break;
case 2: if ((r->infior + INFERIORITY) > wx) {
             // A non-unique keys' inferiority must be gt neighbouring
             // unique keys. Strand bit included, but that shouldn't matter.
             wx ^= r->infior + INFERIORITY;
             kc->kctndx[ndx] ^= wx;
        }
        _incr_at_ndx(kc, r->left, ndx)
        break;
case 3: EPQ(dbg > 4, "2nd+ unique, extend range at %u, infior %u", b2pos, r->infior >> INFIOR_SHFT);
        // A uniqs' inferior must be ge neighbouring uniqs' infior.
        wx ^= r->infior;
        kc->kctndx[ndx] ^= wx; // XXX: why no -ge test before infior replacement?
        ++r->uqct;
        r->left = decr_excise(kc, r->left);
        ASSERT(r->left >= 0, return ~0ul, "[%lu] left? %d", b2pos, r->left);
        if (r->infior <= wx) {
            // special case 4 is for initialization
case 4:     r->infior = wx & R_INFIOR; // only set strand bit below
        }
        r->left = kc->ext;
        _incr_at_ndx(kc, r->left, ndx);
    }
    ++Next_nts->count;
    uint64_t offs = ndx ? kc->kct[ndx-1] : 0;

    // ts.idx of this key minus ts.idx of former is this keys' next-Nts-length.
    r->ct = (kc->kct[ndx] - offs - Next_nts->excise_ct);
    ASSERT(r->ct != 0, return ~0ul);
    r->ct = r->ct == 1ul; // if only one, the key has become uniq

    return offs + add; // return offset to next_nt.
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
    running r = {0};
    r.ct = 4;                   // treat last boundary as a first unique

    EPQ0(dbg > 3, "----[\t%s%s:%u..%u(-%u)\t]----\tdna:", strlen(hdr) < 8 ? "\t":"", 
            hdr, b2pos, next->s, next->s + next->l);
    print_dna(dna, dbg > 3);
    if (dbg > 6) show_list(kc, h->bnd);

    //dbg = 7 * !strncmp(hdr, "GL000207.1", strlen(hdr));
    for(;b2pos < next->s;++b2pos) { // until next event
        // get offset to next Nt
        uint64_t offs = eval_ndx(kc, &r, inter, dna, rc, b2pos);
        ASSERT(offs != ~0ul, return -EFAULT);

        EPQ(dbg > 6, "offs:%lu\tnext Nts (%lu part, %ux, byte %lu(%x)):",
            offs, 4 - (offs & 3), r.ct, offs >> 2, kc->ts[offs >> 2]);
        print_2dna(kc->ts[offs >> 2], kc->ts[offs >> 2] >> ((offs & 3) << 1), dbg > 6);

        uint8_t b2 = (kc->ts[offs >> 2] >> ((offs & 3) << 1)) & 3;
        _verify_seq(b2pos, h->s_s, kc, "[%u]: sb2:%c, got %c, kctndx 0x%lx",
                return print_2dna(dna, rc), b2pos,b6(sb2<<1), b6(b2<<1),
                kc->kctndx[_get_ndx_and_strand(offs, offs, dna, rc)])

        dna = _seq_next(b2, dna, rc);

        EPQ0(dbg > 5, "%u:\t", b2pos); print_2dna(dna, rc, dbg > 5);

        if (r.left) { // within uniq range: keep filling buffer
            if (--r.left) {
                if (r.left == kc->ext - 1) {
                    EPQ(dbg > 6, "setting uniq jump destination");
                    // Sucessively overwritten until region completed - left became 0.
                    inter->dna = dna;
                    inter->l = b2pos - inter->s;
                }
                r.ct |= 2;
            } else { //buf just became full, without extensions.
                EPQ(dbg > 3 && inter->l, "[%u]\t%u - %u\t", kc->bd_l, inter->s, inter->s + inter->l);
                if (inter->l && inter != last) {
                    inter->corr = last->corr;
                    _buf_grow0(kc->bd, 2ul);
                    h->bnd.insert(kc->bdit, kc->bd_l++);
                    last = kc->bd + lastx;
                    next = kc->bd + *kc->bdit;
                }
                inter = kc->bd + kc->bd_l; // no longer last at least
                inter->l = 0;
                r.infior = 0;

                for (unsigned i = kc->ext; i--> 0 ;) {
                    if (update_wlkr(kc, i) < 0)
                        return -EFAULT;
                }
            }
        }
    }
    ASSERT(b2pos == next->s, return -EFAULT, "%u, %u", b2pos, next->s);
    if (r.left) {
        EPQ(dbg > 3, "[%u]\t%u - %u...\t", kc->bd_l, inter->s, next->s + next->l);
        print_2dna(last->dna, inter->at_dna, dbg > 1);

        r.left = decr_excise(kc, r.left);
        ASSERT(r.left >= 0, return -EFAULT, "[%lu] left? %d", b2pos, r.left);
        // if last and next have now become adjoining, merge
        EPQ(dbg > 3, "joining %u(-%u) til %u\t", inter->s,
            inter->s + inter->l, next->s + next->l);
        next->at_dna = inter->at_dna;
        next->l += next->s - inter->s;
        next->s = inter->s;
        next->corr += inter->corr;
 
    }
    return r.uqct;
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
    if (dbg > 5)
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
    for (Walker *w = kc->wlkr; w != kc->wlkr + kc->kct_l; ++w)
        w->count = 0u;

    return uqct;
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


