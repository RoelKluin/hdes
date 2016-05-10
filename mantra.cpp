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
show_mantras(kct_t C*C kc)
{
    unsigned j = 0;

    EPR("-- Mantras --");
    std::list<Mantra>::iterator it = kc->bnd->begin();
    if (it != kc->bnd->end()) {
        do {
            EPR("[%u]:\t(" Pfmt "+) %u - %u", j++, (*it).corr, (*it).s, (*it).e);
        } while (++it != kc->bnd->end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

