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

int
save_boundaries(struct gzfh_t* fhout, Key_t* kc)
{
    int res = -EFAULT;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;

    // 1: buffer sizes
    __WRITE_VAL(kc->h_l, fhout)
    __WRITE_VAL(kc->id_l, fhout)
    __WRITE_VAL(kc->bnd_l, fhout)
    __WRITE_VAL(kc->bnd_m, fhout)

    // 2: next contig id's, because they are somewhat readable.
    __WRITE_PTR(kc->id, fhout, kc->id_l)
    __WRITE_PTR(kc->bnd, fhout, kc->bnd_l)

    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h)
    {
        __WRITE_VAL(h->p_l, fhout)
        __WRITE_VAL(h->len, fhout)
        // sequence, if read, is not stored.
        __WRITE_VAL(h->ido, fhout)
    }
    res = 0;
err:
    rclose(fhout);
    return res;
}

int
load_boundaries(struct gzfh_t* fhin, Key_t* kc)
{
    int res;
    kc->id = NULL;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    // 1: buffer sizes
    __READ_VAL(kc->h_l, fhin) // header size: how many contigs (and uniq regions, fhin)
    __READ_VAL(kc->id_l, fhin)
    __READ_VAL(kc->bnd_l, fhin)
    __READ_VAL(kc->bnd_m, fhin)

    // 2: next contig id's.
    __READ_PTR(kc->id, fhin, kc->id_l)

    kc->bnd = (Mantra*)malloc((1ul << kc->bnd_m) * sizeof(Mantra));
    ASSERT(kc->bnd != NULL, return -ENOMEM);
    if (fhin->read(fhin, (char*)kc->bnd, kc->bnd_l * sizeof(Mantra)))
        return -EFAULT;

    kc->h = (Hdr*)malloc(kc->h_l * sizeof(Hdr));
    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h) {
        __READ_VAL(h->p_l, fhin)
        __READ_VAL(h->len, fhin)
        __READ_VAL(h->ido, fhin)
    }
    res = 0;
err:
    return res;
}

int
save_seqb2(struct gzfh_t* fhout, Key_t* kc)
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
load_seqb2(struct gzfh_t* fhin, Key_t* kc)
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

int save_kc(struct gzfh_t* fhout, Key_t* kc)
{
    uint32_t* buf = NULL;
    int res = -EFAULT;
    ASSERT(fhout->fp == NULL, goto err);
    _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);
    res = -EFAULT;
    __WRITE_VAL(kc->hkoffs_l, fhout)
    __WRITE_VAL(kc->kct_l, fhout)
    __WRITE_VAL(kc->hkoffs_m, fhout)
    __WRITE_PTR(kc->hkoffs, fhout, kc->hkoffs_l)

    __WRITE_PTR(kc->kct, fhout, kc->kct_l)
    res = rclose(fhout);
err:
    if (buf)
        free(buf);
    rclose(fhout);
    return res;
}

int load_kc(struct gzfh_t* fhin, Key_t* kc)
{
    int res;
    _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);
    res = -EFAULT;
    kc->kct = NULL;
    // 0: version number
    // 1: buffer sizes
    __READ_VAL(kc->hkoffs_l, fhin)
    __READ_VAL(kc->kct_l, fhin)
    __READ_VAL(kc->hkoffs_m, fhin)
    kc->hkoffs = (uint32_t*)malloc((1ul << kc->hkoffs_m) * sizeof(uint32_t));
    ASSERT(kc->hkoffs != NULL, return -ENOMEM);
    if (fhin->read(fhin, (char*)kc->hkoffs, kc->hkoffs_l * sizeof(uint32_t)))
        return -EFAULT;

    __READ_PTR(kc->kct, fhin, kc->kct_l);

    for (uint32_t *x = kc->contxt_idx; x != kc->contxt_idx + KEYNT_BUFSZ; ++x)
        *x = kc->kct_l;

    uint8_t *s = kc->s;
    Hdr *h = kc->h;
    uint32_t *hkoffs = kc->hkoffs;
    for (uint32_t *k = kc->kct; k != kc->kct + kc->kct_l; ++k) {

        while (k >= kc->kct + *hkoffs) {
            EPR0("loaded kcs:%s", kc->id + h->ido);
            s += ++h->len;
            ++hkoffs;
        }
        // construct index from sequence
        keyseq_t seq = {.p = *k};
        uint32_t ndx = build_ndx_kct(kc, seq, s);
        NB(ndx < KEYNT_BUFSZ);

        kc->contxt_idx[ndx] = k - kc->kct;
    }
    res = 0;
err:
    return res;
}

int ammend_kc(struct gzfh_t* fhin, Key_t* kc)
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

