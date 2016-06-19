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
            while (h - kc->h != (*it).ho)
                ++h;

            EPR("[%u%c]:\t>%s (%u+)%u - %u", j++, it==here?'*':' ',
                    kc->id + h->ido, (*it).corr, (*it).s >>1, (*it).e >> 1);
        } while (++it != kc->bnd->end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

