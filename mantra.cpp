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

#include <list>
#include "fa.h"

int
show_mantras(kct_t C*C kc, std::list<Mantra>::iterator here)
{
    unsigned j = 0;

    EPR("-- Mantras --");
    Hdr* h = kc->h;
    std::list<Mantra>::iterator it = kc->bnd->begin();
    if (it != kc->bnd->end()) {
        do {
            uint32_t ke = (*it).ke;
            if (ke <= kc->hkoffs[kc->h_l-1]) {
                while (ke > kc->hkoffs[h - kc->h])
                    ++h;
            } else {
                h = kc->h;
                for (uint32_t* hkoffs = kc->hkoffs + kc->h_l;
                        hkoffs != (kc->hkoffs + kc->hkoffs_l) && ke >= *hkoffs; ++hkoffs, ++h)
                    if (h == kc->h + kc->h_l)
                        h = kc->h;
            }
            //EPR("%s",  ke != kc->hkoffs[h - kc->h] ? "kepos(kc, it) - 2" : "h->end");
            uint32_t end = ke != kc->hkoffs[h - kc->h] ? kepos(kc, it) - 2 : h->end;

            EPR("[%u%c]:\t>%s (%u+)%u - %u\t(ke:%u)", j++, it==here?'*':' ',
                    kc->id + h->ido, (*it).corr, (*it).s >>1, end >> 1, ke);
            ++h;
        } while (++it != kc->bnd->end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

