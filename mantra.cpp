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
            uint32_t end = ke != kc->hkoffs[h - kc->h] ? kepos(kc, it) - 2 : h->end;

            EPR("[%u%c]:\t>%s (%u+)%u - %u\t(ke:%u)", j++, it==here?'*':' ',
                    kc->id + h->ido, (*it).corr, (*it).s >>1, end >> 1, ke);
            if (ke >= kc->hkoffs[h - kc->h])
                ++h;
        } while (++it != kc->bnd->end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

