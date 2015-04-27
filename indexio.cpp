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


int write1(struct gzfh_t* fhout, kct_t* kc)
{
    int ret = -EFAULT;
    if (fhout->write(fhout, (const char*)&kc->bd_l, sizeof(uint32_t)) < 0)
        return ret;
    EPR("%u", kc->bd_l);
    if (fhout->write(fhout, (const char*)&kc->id_l, sizeof(uint32_t)) < 0)
        return ret;
    if (fhout->write(fhout, (const char*)kc->id, kc->id_l * sizeof(char)) < 0)
        return ret;
    if (fhout->write(fhout, (const char*)&kc->kct_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t len = kc->h.size();
    if (fhout->write(fhout, (const char*)&len, sizeof(uint32_t)) < 0)
        return ret;
    for(std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h)
    {
        EPR("%s..", kc->id + (*h)->part[0]);
        if (fhout->write(fhout, (const char*)&(*h)->end_pos, sizeof(uint32_t)) < 0)
            return ret;
        uint32_t len = (*h)->bnd.size();
        if (fhout->write(fhout, (const char*)&len, sizeof(uint32_t)) < 0)
            return ret;
        for(std::list<uint32_t>::iterator b = (*h)->bnd.begin(); b != (*h)->bnd.end(); ++b)
        {
            if (fhout->write(fhout, (const char*)&*b, sizeof(uint32_t)) < 0)
                return ret;
        }
        // sequence is not stored.
        if (fhout->write(fhout, (const char*)(*h)->part, 10 * sizeof(uint16_t)) < 0)
            return ret;
        if (fhout->write(fhout, (const char*)&(*h)->p_l, sizeof(uint8_t)) < 0)
            return ret;
        // p_l is not needed anymore
    }
    // FIXME: full flush silences valgrind errors on close, about uninit
    // values - either when writing buf or kc->kct.
    ASSERT(gzflush(fhout->io, Z_FULL_FLUSH) == Z_OK, return -EFAULT);
    EPR("after flush1");
    if (fhout->write(fhout, (const char*)kc->kct, kc->kct_l * sizeof(Kct)) < 0)
        return ret;

    len = kc->kct_l * 2;
    uint32_t* buf = (uint32_t*)malloc(sizeof(uint32_t) * len);
    if (buf == NULL) return -ENOMEM;
    uint32_t i = 0;
    // FIXME: may need 64 bit index.
    for (uint32_t ndx = 0u; i != len; ++ndx) {
        if (kc->kcsndx[ndx] >= kc->kct_l) {
            kc->kcsndx[ndx] = kc->kct_l; // make high bit range available
        } else {
            buf[i++] = ndx;
            buf[i++] = kc->kcsndx[ndx];
        }
    }
    if (fhout->write(fhout, (const char*)buf, sizeof(uint32_t) * len) < 0)
        goto err;
    ASSERT(gzflush(fhout->io, Z_FULL_FLUSH) == Z_OK, return -EFAULT);
    EPR("after flush2");

    len = kc->bd_l * sizeof(Bnd);
    // when kc->bd is not packed, problem at 0th.
    /*for (unsigned i = 0; i != kc->bd_l; ++i) {
        Bnd* bd = &kc->bd[i];
        EPR("%u:\t%u\t%u\t%u\t[%lu]", i, bd->s, bd->l, bd->i, sizeof(kc->bd[i]));
        print_dnarc(bd->at_dna, bd->dna);
    }*/
    // kc->wlkr and kc->wbuf aren't filled yet.
    if (fhout->write(fhout, (const char*)kc->bd, len) < 0)
        return ret;
    ASSERT(gzflush(fhout->io, Z_FULL_FLUSH) == Z_OK, return -EFAULT);
    EPR("after flush3");
    ret = 0;
err:
    free(buf);
    return ret;
}

int restore1(struct gzfh_t* fhin, kct_t* kc)
{
    int ret = -EFAULT;
    kc->bd = NULL; kc->id = NULL; kc->kct = NULL;
    if (fhin->read(fhin, (char*)&kc->bd_l, sizeof(uint32_t)) < 0)
        return ret;
    EPR("%u", kc->bd_l);
    if (fhin->read(fhin, (char*)&kc->id_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t len = kc->id_l * sizeof(char);
    kc->id = (char*)malloc(len);
    if (fhin->read(fhin, (char*)kc->id, len) < 0)
        return ret;
    if (fhin->read(fhin, (char*)&kc->kct_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t blen, val;
    if (fhin->read(fhin, (char*)&len, sizeof(uint32_t)) < 0)
        return ret;
    for (uint32_t i=0; i != len; ++i) {
        Hdr *h = new Hdr;
        h->s = NULL; h->s_l = 0u; // XXX: may be removed in future.
        if (fhin->read(fhin, (char*)&h->end_pos, sizeof(uint32_t)) < 0)
            return ret;
        if (fhin->read(fhin, (char*)&blen, sizeof(uint32_t)) < 0)
            return ret;
        for (uint32_t j=0; j != blen; ++j) {
            if (fhin->read(fhin, (char*)&val, sizeof(uint32_t)) < 0)
                return ret;
            h->bnd.push_back(val);
        }
        if (fhin->read(fhin, (char*)h->part, 10 * sizeof(uint16_t)) < 0)
            return ret;
        if (fhin->read(fhin, (char*)&h->p_l, sizeof(uint8_t)) < 0)
            return ret;
        kc->h.push_back(h);
        kc->hdr.insert(std::pair<char*, Hdr*>(kc->id + h->part[ID], kc->h.back()));
    }
    len = kc->kct_l * sizeof(Kct);
    kc->kct = (Kct*)malloc(len);
    if (fhin->read(fhin, (char*)kc->kct, len) < 0)
        return ret;
    memset(kc->kcsndx, kc->kct_l, KEYNT_BUFSZ * sizeof(kc->kcsndx[0]));
    for (uint32_t i=0; i != kc->kct_l; ++i) {
        if (fhin->read(fhin, (char*)&val, sizeof(uint32_t)) < 0)
            return ret;
        ASSERT(val < KEYNT_BUFSZ, return -EFAULT, "%u/%u: %u > KEYNT_BUFSZ", i, kc->kct_l, val);
        if (fhin->read(fhin, (char*)&len, sizeof(uint32_t)) < 0)
            return ret;
        kc->kcsndx[val] = len;
    }
    len = next_pow2(kc->bd_l + 2);
    kc->bd_m = __builtin_clz(len);
    kc->bd = (Bnd*)malloc(len);
    len = kc->bd_l * sizeof(Bnd);
    if (fhin->read(fhin, (char*)kc->bd, len) < 0)
        return ret;
    return 0;
}

