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

int
show_mantras(key_t C*C kc, Mantra* obnd, unsigned obnd_l, Mantra* at)
{
    unsigned j = 0;

    EPR("-- Mantras --");
    if (obnd_l) {
        Hdr* h = kc->h;
        Mantra* bnd = kc->bnd_l ? kc->bnd : at;
        while (bnd != obnd + obnd_l) {

            NB (h - kc->h <= bnd->ho);
            while (h - kc->h != bnd->ho)
                ++h;

            EPR("[%u%c]:\t>%s (%u+)%u - %u", j++, bnd==at?'*':' ',
                    kc->id + h->ido, bnd->corr, bnd->s >>1, bnd->e >> 1);

            if (++bnd == kc->bnd + kc->bnd_l) {
                bnd = at;
                if (bnd - obnd == obnd_l)
                    break;
            }
        }
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

