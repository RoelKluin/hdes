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
show_mantras(kct_t C*C kc, Hdr *C h)
{
    unsigned j = 0;

    std::list<Mantra>::iterator it = h->bnd.begin();
    if (it != h->bnd.end()) {
        do {
            keyseq_t seq = {0};
            pos_t p = (*it).s;
            seq.p = p + KEY_WIDTH;
            _build_key(kc, seq, p, seq.p, seq.t);
            EPR0("[%s%u]:\t(" Pfmt "+) %lu - %lu\t", it != kc->bdit ? "" : "* ",
                    j++, (*it).corr, (*it).s + h->s_s, (*it).e + h->s_s);
            print_dna(seq.dna);
        } while (++it != h->bnd.end());
    } else {
        EPR("[ entirely mapable ]");
    }
    return -1;
}

