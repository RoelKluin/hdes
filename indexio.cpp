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
#include "fa.h"
#include "gz.h"

// XXX: these need to be updated
// read write of 64 bits may require some more work.

#define __WRITE_VAL(x) \
    ASSERT(fhout->fp != NULL, return -EFAULT);\
    ASSERT(fhout->write != NULL, return -EFAULT);\
    EPQ(dbg > 2, "writing value for " #x "");\
    if (fhout->write(fhout, (const char*)&(x), sizeof(x)) < 0)\
        goto err;
#define __WRITE_PTR(x, l) \
    EPQ(dbg > 2, "writing poiner " #x "");\
    if (fhout->write(fhout, (const char*)(x), (l) * sizeof(*(x))) < 0)\
        goto err;

#define __READ_VAL(x) \
    ASSERT(fhin->fp != NULL, return -EFAULT);\
    ASSERT(fhin->read != NULL, return -EFAULT);\
    if (fhin->read(fhin, (char*)&(x), sizeof(x)) < 0) \
        return res;
#define __READ_PTR(x, l) \
    x = (decltype(x))malloc((l) * sizeof(*(x)));\
    ASSERT(x != NULL, return -ENOMEM);\
    if (fhin->read(fhin, (char*)x, l * sizeof(*(x))))\
        return res;

int
save_boundaries(struct gzfh_t* fhout, kct_t* kc)
{
    int res = -EFAULT;
    uint32_t val, len = kc->h.size();
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;
    // 1: buffer sizes
    __WRITE_VAL(len)
    __WRITE_VAL(kc->id_l)
    __WRITE_VAL(kc->bd_l)

    // 2: next contig id's, because they are somewhat readable.
    __WRITE_PTR(kc->id, kc->id_l)
    __WRITE_PTR(kc->bd, kc->bd_l)

    for(std::list<Hdr*>::iterator hit = kc->h.begin(); hit != kc->h.end(); ++hit)
    {
        Hdr* h = *hit;
        __WRITE_VAL(h->end_pos)
        len = h->bnd.size();
        __WRITE_VAL(len)
        __WRITE_VAL(h->s_s)
        for(std::list<uint32_t>::iterator b = h->bnd.begin(); b != h->bnd.end(); ++b) {
            val = *b; // NB: calls list's overloaded function: get stored value.
            __WRITE_VAL(val)
        }
        __WRITE_VAL(h->p_l)
        // sequence, if read, is not stored.
        __WRITE_PTR(h->part, h->p_l)
    }
    res = 0;
err:
    rclose(fhout);
    return res;
}

int
load_boundaries(struct gzfh_t* fhin, kct_t* kc)
{
    int res;
    uint32_t len, blen, val;
    kc->bd = NULL;
    kc->id = NULL;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    // 1: buffer sizes
    __READ_VAL(len) // header size: how many contigs (and uniq regions)
    __READ_VAL(kc->id_l)
    __READ_VAL(kc->bd_l)

    // 2: next contig id's.
    __READ_PTR(kc->id, kc->id_l)

    // special: we reserve more space than stored so we can continue using it.
    kc->bd_m = __builtin_ctz(next_pow2(kc->bd_l + 2));
    blen = kc->bd_l * sizeof(Bnd);
    kc->bd = (Bnd*)malloc((1 << kc->bd_m) * sizeof(Bnd));
    ASSERT(kc->bd != NULL, return -ENOMEM);
    if (fhin->read(fhin, (char*)kc->bd, blen) < 0)
        return res;

    for (uint32_t i=0; i != len; ++i) {
        Hdr *h = new Hdr;
        __READ_VAL(h->end_pos)
        __READ_VAL(blen)
        __READ_VAL(h->s_s)
        for (uint32_t j=0; j != blen; ++j) {
            __READ_VAL(val)
            h->bnd.push_back(val);
        }
        __READ_VAL(h->p_l)
        __READ_PTR(h->part, h->p_l)
        kc->h.push_back(h);
    }
    res = 0;
err:
    return res;
}

int
save_seqb2(struct gzfh_t* fhout, kct_t* kc)
{
    int res = -EFAULT;
    uint64_t len64;
    ASSERT(fhout->fp == NULL, goto err, "%s", fhout->name);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    __WRITE_VAL(kc->s_l)
    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __WRITE_PTR(kc->s, len64)
    res = 0;
err:
    rclose(fhout);
    return res;
}

int
load_seqb2(struct gzfh_t* fhin, kct_t* kc)
{
    int res;
    uint64_t len64;
    kc->s = NULL;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    __READ_VAL(kc->s_l)
    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __READ_PTR(kc->s, len64)
    res = 0;
err:
    return res;
}

int save_nextnts(struct gzfh_t* fhout, kct_t* kc)
{
    int res = -EFAULT;
    uint64_t len64;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    __WRITE_VAL(kc->ts_l)
    len64 = (kc->ts_l >> 2) + !!(kc->ts_l & 3);
    __WRITE_PTR(kc->ts, len64)
    res = 0;
err:
    rclose(fhout);
    return res;
}
int load_nextnts(struct gzfh_t* fhin, kct_t* kc)
{
    int res = -EFAULT;
    uint64_t len64;
    kc->ts = NULL;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    __READ_VAL(kc->ts_l)
    len64 = (kc->ts_l >> 2) + !!(kc->ts_l & 3);
    __READ_PTR(kc->ts, len64)
    res = 0;
err:
    return res;

}

int save_kc(struct gzfh_t* fhout, kct_t* kc)
{
    uint64_t* buf = NULL;
    int res = -EFAULT;
    uint32_t i;
    uint64_t len64;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    __WRITE_VAL(kc->kct_l)
    __WRITE_PTR(kc->kct, kc->kct_l)
    len64 = kc->kct_l * 2;
    buf = (uint64_t*)malloc(sizeof(uint64_t) * len64);
    ASSERT(buf != NULL, return -ENOMEM);
    i = 0;
    // TODO: put this conversion in kct_convert and only write kc->kctndx here..
    for (uint64_t ndx = 0ul; ndx != KEYNT_BUFSZ; ++ndx) {

        // ..then INDEX_MASK'ing isn't necessary test should be
        // kc->kctndx[ndx] == ~0ul
        if ((kc->kctndx[ndx] & INDEX_MASK) >= kc->kct_l) {
            kc->kctndx[ndx] = kc->kct_l; // make high bit range available
        } else {
            buf[i++] = ndx;
            buf[i++] = kc->kctndx[ndx];
        }
    }
    __WRITE_PTR(buf, len64)
    res = rclose(fhout);
err:
    if (buf)
        free(buf);
    rclose(fhout);
    return res;
}

int load_kc(struct gzfh_t* fhin, kct_t* kc)
{
    int res;
    uint64_t val64, len64;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    kc->kct = NULL;
    // 0: version number
    // 1: buffer sizes
    __READ_VAL(kc->kct_l)
    __READ_PTR(kc->kct, kc->kct_l);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->kctndx[i] = kc->kct_l;
    for (uint64_t i=0ul; i != kc->kct_l; ++i) {
        __READ_VAL(val64)
        ASSERT(val64 < KEYNT_BUFSZ, return -EFAULT, "%lu/%lu: %lu > KEYNT_BUFSZ", i, kc->kct_l, val64);
        __READ_VAL(len64)
        kc->kctndx[val64] = len64;
    }
    res = 0;
err:
    return res;
}

