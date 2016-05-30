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
            EPR("[%u%c]:\t>%s (%u+)%u - %u\t(ke:%u)", j++, it==here?'*':' ',
                    kc->id + h->ido, (*it).corr, (*it).s >>1,
                    ((*it).ke == kc->hkoffs[h - kc->h] ? h->end: kepos(kc, it)) >>1, (*it).ke);
            if ((*it).ke == kc->hkoffs[h - kc->h])
                ++h;
        } while (++it != kc->bnd->end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

