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

#define __WRITE_INIT(fhout) \
    do {\
        ASSERT(fhout, goto err);\
        ASSERT(fhout->fp == NULL, goto err, "%s", fhout->name);\
        _ACTION(set_io_fh(fhout, 1), "opening %s for writing", fhout->name);\
        ASSERT(fhout->fp, goto err, "fp:%s", fhout->name);\
        ASSERT(fhout->write, goto err, "write:%s", fhout->name);\
    } while (0)

#define __READ_INIT(fhin) \
    do {\
        ASSERT(fhin, goto err);\
        _ACTION(set_io_fh(fhin, 2), "opening %s for reading", fhin->name);\
        ASSERT(fhin->fp, goto err, "fp:%s", fhin->name);\
        ASSERT(fhin->read, goto err, "read:%s", fhin->name);\
    } while (0)

#define __WRITE_VAL(x, fhout) \
    do {\
        if (fhout->write(fhout, (const char*)&(x), sizeof(x)) < 0)\
            goto err;\
    } while (0)

#define __WRITE_PTR(x, fhout, l) \
    do {\
        if (fhout->write(fhout, (const char*)(x), (l) * sizeof(*(x))) < 0)\
            goto err;\
    } while (0)

#define __READ_VAL(x, fhin) \
    do {\
        if (fhin->read(fhin, (char*)&(x), sizeof(x)) < 0) \
            return -EFAULT;\
    } while (0)

#define __READ_PTR_(x, fhin, l, m) \
    do {\
        x = (decltype(x))malloc((m) * sizeof(*(x)));\
        ASSERT(x != NULL, return -ENOMEM);\
        if (fhin->read(fhin, (char*)x, l * sizeof(*(x))))\
            return -EFAULT;\
    } while (0)

#define __READ_PTR(x, fhin, l) __READ_PTR_(x, fhin, l, l)
#define __READ_LMPTR(x, fhin) __READ_PTR_(x, fhin, x##_l, (1ul << x##_m))


int
save_seqb2(struct gzfh_t* fhout, Key_t* kc)
{
    int res = -EFAULT;
    uint64_t len64;
    __WRITE_INIT(fhout);

    __WRITE_VAL(kc->h_l, fhout);
    __WRITE_VAL(kc->id_l, fhout);
    __WRITE_VAL(kc->s_l, fhout);

    __WRITE_PTR(kc->id, fhout, kc->id_l);
    __WRITE_PTR(kc->h, fhout, kc->h_l);

    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __WRITE_PTR(kc->s, fhout, len64);
    res = 0;
err:
    rclose(fhout);
    return res;
}

int
load_seqb2(struct gzfh_t* fhin, Key_t* kc)
{
    int res = -EFAULT;
    uint64_t len64;
    __READ_INIT(fhin);

    __READ_VAL(kc->h_l, fhin); // header size: how many contigs (and uniq regions, fhin)
    __READ_VAL(kc->id_l, fhin);
    __READ_VAL(kc->s_l, fhin);

    __READ_PTR(kc->id, fhin, kc->id_l);
    __READ_PTR(kc->h, fhin, kc->h_l);

    len64 = (kc->s_l >> 2) + !!(kc->s_l & 3);
    __READ_PTR(kc->s, fhin, len64);
    res = 0;
err:
    return res;
}

int save_kc(struct gzfh_t* fhout, Key_t* kc)
{
    int res = -EFAULT;
    __WRITE_INIT(fhout);

    __WRITE_VAL(kc->kct_l, fhout);
    __WRITE_VAL(kc->ext_iter_l, fhout);
    __WRITE_VAL(kc->bnd_l, fhout);
    __WRITE_VAL(kc->hkoffs_l, fhout);

    __WRITE_VAL(kc->kct_m, fhout);
    __WRITE_VAL(kc->ext_iter_m, fhout);
    __WRITE_VAL(kc->bnd_m, fhout);
    __WRITE_VAL(kc->hkoffs_m, fhout);

    __WRITE_PTR(kc->kct, fhout, kc->kct_l);
    __WRITE_PTR(kc->ext_iter, fhout, kc->ext_iter_l);
    __WRITE_PTR(kc->bnd, fhout, kc->bnd_l);
    __WRITE_PTR(kc->hkoffs, fhout, kc->hkoffs_l);

    res = rclose(fhout);
err:
    rclose(fhout);
    return res;
}
// TODO: specialization dependent on kc->ext not 0?
int load_kc(struct gzfh_t* fhin, Key_t* kc)
{
    int res = -EFAULT;
    uint8_t *s = kc->s;
    Hdr *h = kc->h;
    uint32_t *hkoffs;
    __READ_INIT(fhin);

    __READ_VAL(kc->kct_l, fhin);
    __READ_VAL(kc->ext_iter_l, fhin);
    __READ_VAL(kc->bnd_l, fhin);
    __READ_VAL(kc->hkoffs_l, fhin);

    __READ_VAL(kc->kct_m, fhin);
    __READ_VAL(kc->ext_iter_m, fhin);
    __READ_VAL(kc->bnd_m, fhin);
    __READ_VAL(kc->hkoffs_m, fhin);

    __READ_LMPTR(kc->kct, fhin);
    __READ_LMPTR(kc->ext_iter, fhin);
    __READ_LMPTR(kc->bnd, fhin);
    __READ_LMPTR(kc->hkoffs, fhin);

    for (uint32_t *x = kc->contxt_idx; x != kc->contxt_idx + KEYNT_BUFSZ; ++x)
        *x = kc->kct_l;

    hkoffs = kc->hkoffs;
    for (uint32_t *k = kc->kct; k != kc->kct + kc->kct_l; ++k) {

        while (K_OFFS(kc, k) == *hkoffs) {
            if (h - kc->h != kc->h_l - 1) {
                s += h++->len;
            } else {
                h = kc->h;
                s = kc->s;
            }
            ++hkoffs;
        }
        // construct index from sequence
        // there can be keys with a dupbit set. They have a position, but that's just one of
        // many possible.
        keyseq_t seq = {.p = *k & ~DUP_BIT};
        NB((seq.p >> 1) <= h->end, "%x", seq.p);
        uint32_t ndx = build_ndx_kct(kc, seq, s);
        NB(ndx < KEYNT_BUFSZ);

        kc->contxt_idx[ndx] = K_OFFS(kc, k);
    }
    res = 0;
err:
    return res;
}

