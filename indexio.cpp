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
    EPQ(dbg > 5, "writing value for " #x ":0x%lx", (uint64_t)(x));\
    if (fhout->write(fhout, (const char*)&(x), sizeof(x)) < 0)\
        goto err;
#define __WRITE_PTR(x, l) \
    EPQ(dbg > 5, "writing poiner " #x "");\
    if (fhout->write(fhout, (const char*)(x), (l) * sizeof(*(x))) < 0)\
        goto err;
#define __WRITE_LMPTR(buf)\
    __WRITE_VAL(buf##_l)\
    __WRITE_VAL(buf##_m)\
    __WRITE_PTR(buf, buf##_l)

#define __READ_VAL(x) \
    ASSERT(fhin->fp != NULL, return -EFAULT);\
    ASSERT(fhin->read != NULL, return -EFAULT);\
    if (fhin->read(fhin, (char*)&(x), sizeof(x)) < 0) \
        return -EFAULT;
#define __READ_PTR(x, l) \
    x = (decltype(x))malloc((l) * sizeof(*(x)));\
    ASSERT(x != NULL, return -ENOMEM);\
    if (fhin->read(fhin, (char*)x, l * sizeof(*(x))))\
        return -EFAULT;
#define __READ_LMPTR(buf) \
    __READ_VAL(buf##_l)\
    __READ_VAL(buf##_m)\
    buf = (decltype(buf))malloc((1ul << buf##_m) * sizeof(*(buf)));\
    ASSERT(buf != NULL, return -ENOMEM);\
    if (fhin->read(fhin, (char*)buf, buf##_l * sizeof(*(buf))))\
        return -EFAULT;

int
save_boundaries(struct gzfh_t* fhout, kct_t* kc)
{
    int res = -EFAULT;
    uint32_t val, len;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    // 1: buffer sizes
    __WRITE_VAL(kc->h_l)
    __WRITE_VAL(kc->id_l)

    // 2: next contig id's, because they are somewhat readable.
    __WRITE_PTR(kc->id, kc->id_l)

    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h)
    {
        __WRITE_VAL(h->end_pos)
        len = h->bnd->size();
        __WRITE_VAL(len)
        __WRITE_VAL(h->s_s)
        for(std::list<Mantra>::iterator b = h->bnd->begin(); b != h->bnd->end(); ++b) {
            val = (*b).corr;
            __WRITE_VAL(val)
            val = (*b).s;
            __WRITE_VAL(val)
            val = (*b).e;
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
    uint32_t blen;
    kc->id = NULL;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    // 1: buffer sizes
    __READ_VAL(kc->h_l) // header size: how many contigs (and uniq regions)
    __READ_VAL(kc->id_l)

    // 2: next contig id's.
    __READ_PTR(kc->id, kc->id_l)
    kc->h = (Hdr*)malloc(kc->h_l * sizeof(Hdr));

    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h) {
        __READ_VAL(h->end_pos)
        __READ_VAL(blen)
        __READ_VAL(h->s_s)
        for (uint32_t j=0; j != blen; ++j) {
            Mantra contig = {0};
            __READ_VAL(contig.corr)
            __READ_VAL(contig.s)
            __READ_VAL(contig.e)
            h->bnd = new std::list<Mantra>();
            h->bnd->push_back(contig);
        }
        __READ_VAL(h->p_l)
        __READ_PTR(h->part, h->p_l)
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

int save_kc(struct gzfh_t* fhout, kct_t* kc)
{
    uint32_t* buf = NULL;
    int res = -EFAULT;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;
    __WRITE_LMPTR(kc->hk)

    __WRITE_VAL(kc->kct_l)
    __WRITE_PTR(kc->kct, kc->kct_l)
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
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    kc->kct = NULL;
    // 0: version number
    // 1: buffer sizes
    __READ_LMPTR(kc->hk)

    __READ_VAL(kc->kct_l)
    __READ_PTR(kc->kct, kc->kct_l);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->ndxkct[i] = kc->kct_l;
    for (unsigned i=0u; i != kc->kct_l; ++i) {
        //ASSERT(kc->kct[i] < KEYNT_BUFSZ, return -EFAULT, "%u/%u: %lu > KEYNT_BUFSZ(%lu)",
        //        i, kc->kct_l, kc->kct[i], KEYNT_BUFSZ);
        kc->ndxkct[kc->kct[i]] = i;
    }
    res = 0;
err:
    return res;
}

int ammend_kc(struct gzfh_t* fhin, kct_t* kc)
{
    int res;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    // 0: version number
    // 1: buffer sizes
    __READ_VAL(kc->kct_l)
    if (fhin->read(fhin, (char*)kc->kct, kc->kct_l * sizeof(*(kc->kct))) < 0)
        return -EFAULT;

    res = 0;
err:
    return res;
}

