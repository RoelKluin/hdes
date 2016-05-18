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


#define __WRITE_VAL(x, fhout) \
    ASSERT(fhout->fp != NULL, return -EFAULT);\
    ASSERT(fhout->write != NULL, return -EFAULT);\
    if (fhout->write(fhout, (const char*)&(x), sizeof(x)) < 0)\
        goto err;
#define __WRITE_PTR(x, fhout, l) \
    if (fhout->write(fhout, (const char*)(x), (l) * sizeof(*(x))) < 0)\
        goto err;
#define __WRITE_LMPTR(buf, fhout)\
    __WRITE_VAL(buf##_l, fhout)\
    __WRITE_VAL(buf##_m, fhout)\
    __WRITE_PTR(buf, fhout, buf##_l)

#define __READ_VAL(x, fhin) \
    ASSERT(fhin->fp != NULL, return -EFAULT);\
    ASSERT(fhin->read != NULL, return -EFAULT);\
    if (fhin->read(fhin, (char*)&(x), sizeof(x)) < 0) \
        return -EFAULT;
#define __READ_PTR(x, fhin, l) \
    x = (decltype(x))malloc((l) * sizeof(*(x)));\
    ASSERT(x != NULL, return -ENOMEM);\
    if (fhin->read(fhin, (char*)x, l * sizeof(*(x))))\
        return -EFAULT;
#define __READ_LMPTR(buf, fhin) \
    __READ_VAL(buf##_l, fhin)\
    __READ_VAL(buf##_m, fhin)\
    buf = (decltype(buf))malloc((1ul << buf##_m) * sizeof(*(buf)));\
    ASSERT(buf != NULL, return -ENOMEM);\
    if (fhin->read(fhin, (char*)buf, buf##_l * sizeof(*(buf))))\
        return -EFAULT;

int
save_boundaries(struct gzfh_t* fhout, kct_t* kc)
{
    int res = -EFAULT;
    uint32_t val;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    // 1: buffer sizes
    __WRITE_VAL(kc->h_l, fhout)
    __WRITE_VAL(kc->id_l, fhout)
    val = kc->bnd->size();
    __WRITE_VAL(val, fhout)

    // 2: next contig id's, because they are somewhat readable.
    __WRITE_PTR(kc->id, fhout, kc->id_l)

    for(std::list<Mantra>::iterator b = kc->bnd->begin(); b != kc->bnd->end(); ++b) {
        val = (*b).corr;
        __WRITE_VAL(val, fhout)
        val = (*b).s;
        __WRITE_VAL(val, fhout)
        val = (*b).e;
        __WRITE_VAL(val, fhout)
        val = (*b).ke;
        __WRITE_VAL(val, fhout)
    }

    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h)
    {
        __WRITE_VAL(h->p_l, fhout)
        // sequence, if read, is not stored.
        __WRITE_VAL(h->ido, fhout)
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
    __READ_VAL(kc->h_l, fhin) // header size: how many contigs (and uniq regions, fhin)
    __READ_VAL(kc->id_l, fhin)
    __READ_VAL(blen, fhin)

    // 2: next contig id's.
    __READ_PTR(kc->id, fhin, kc->id_l)
    kc->bnd = new std::list<Mantra>();

    while (blen--) {
        Mantra m = {0};
        __READ_VAL(m.corr, fhin)
        __READ_VAL(m.s, fhin)
        __READ_VAL(m.e, fhin)
        __READ_VAL(m.ke, fhin)
        kc->bnd->push_front(m);
    }

    kc->h = (Hdr*)malloc(kc->h_l * sizeof(Hdr));
    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h) {
        __READ_VAL(h->p_l, fhin)
        __READ_VAL(h->ido, fhin)
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

    __WRITE_VAL(kc->s_l, fhout)
    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __WRITE_PTR(kc->s, fhout, len64)
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
    __READ_VAL(kc->s_l, fhin)
    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __READ_PTR(kc->s, fhin, len64)
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
    __WRITE_LMPTR(kc->hk, fhout)

    __WRITE_VAL(kc->kct_l, fhout)
    __WRITE_PTR(kc->kct, fhout, kc->kct_l)
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
    __READ_LMPTR(kc->hk, fhin)

    __READ_VAL(kc->kct_l, fhin)
    __READ_PTR(kc->kct, fhin, kc->kct_l);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->contxt_idx[i] = kc->kct_l;

    HK *hk = kc->hk;
    uint8_t *s = kc->s;

    for (uint64_t i=0ul; i != kc->kct_l; ++i) {

        uint32_t *k = kc->kct + i;

        // XXX: ensure the k's contigs are sorted or this won't be efficient.
        while (k >= kc->kct + hk->koffs) {
            s += hk->len;
            NB(hk != kc->hk + kc->hk_l);
            ++hk;
        }
        // construct index from sequence
        keyseq_t seq = {.p = b2pos_of(kc, k)};
        NB(seq.p > NT_WIDTH - 2);

        uint32_t ndx = build_ndx_kct(kc, seq, s);
        NB(ndx < KEYNT_BUFSZ);

        kc->contxt_idx[ndx] = i;
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
    __READ_VAL(kc->kct_l, fhin)
    if (fhin->read(fhin, (char*)kc->kct, kc->kct_l * sizeof(*(kc->kct))) < 0)
        return -EFAULT;

    res = 0;
err:
    return res;
}

