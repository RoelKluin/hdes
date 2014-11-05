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
#include <stdlib.h> // realloc()
#include <ctype.h> // isspace()
#include <errno.h> // ENOMEM
#include <string.h> // memset()
#include "fa.h"

int fa_ndx2(seqb2_t *fa)
{
    uint8_t *s = fa->s;
    uint32_t *const lookup = fa->lookup;
    const unsigned readlength = fa->readlength;
    unsigned key_ct = fa->key_ct;
    const unsigned fa_ent_max = SEQ_MAX_NAME_ETC;
    void* g;
    int (*gc) (void*);
    if (fa->fh[2].io) {
        g = fa->fh[2].io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fa->fh[2].fp;
        gc = (int (*)(void*))&fgetc;
    }
    register uint8_t *b = s + fa->s_l;
    register int c;
    const unsigned width = readlength - KEY_WIDTH;
    uint64_t rot[width];

    while ((c = gc(g)) != '>')
        if (c == '@' && c == -1) goto out; // skip to first
    do {
        assert((uint64_t)fa->s_l + fa_ent_max < (1ul << 32));
        b = s + fa->s_l;
        uint8_t* h = b;
        while (!isspace(c = gc(g)) && (c >= 0)) *b++ = c; // header
        while (c != '\n' && (c >= 0)) c = gc(g);          // comment
        *b = '\0';
        fprintf(stderr, "%s...\n", (char*)h);
        do { // repeat foreach stretch of N's
            register unsigned i;
            uint64_t pos = h - s; // store location of reference string
            for (i = BUF_OFFSET_BYTES; i != 0; --i, pos >>= 8) *++b = pos & 0xff;

            while(isspace(c = gc(g)) || (++i, c == 'N')) {}
            if (i) fprintf(stderr, "Skipped %d N's on %s\n", i, h);

            pos = i; // count of initial Ns
            for (i = NS_OFFSET_BYTES; i != 0; --i, pos >>= 8) *++b = pos & 0xff;

            // (re)initialize key
            register uint64_t dna = 0ul, rev = 0ul;
            do {
                c = b6(c);
                c = -isb6(c) & (c >> 1); // zero [min] for N - but also for A.

                dna = (dna << 2) | c;
                rev = ((uint64_t)c << KEYNT_TOP) | (rev >> 2); // xor 0x2 after loop
                while(isspace(c = gc(g))) {} // get next Nt

                // seq too short for mapping => skip
                if_ever (c == '>' || c == -1) {
                    b -= BUF_OFFSET_BYTES + NS_OFFSET_BYTES; // remove offsets
                    break;
                }
                if ((i & 3) == 0) *++b = '\0';
                *b |= c << ((i & 3) << 1);
            } while (++i != KEY_WIDTH);

            pos += KEY_WIDTH;
            if_ever (c == '>' || c == -1) break;
            rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement

            memset(rot, '\0', readlength * sizeof(uint64_t));
            // seq or revcmp according to 2nd bit of central Nt
            uint64_t lmax = (dna & KEYNT_STRAND) ? rev : dna;
            rot[0] = lmax | (dna & KEYNT_STRAND);
            uint64_t mx = lmax ^ (lmax >> 2);
            unsigned maxi = i = 0;

            lmax = ((lmax >> 1) & ~HALF_KEYNT_MASK) | (lmax & HALF_KEYNT_MASK); // excise out 2nd cNt bit
            lmax &= KEYNT_TRUNC_MASK; // truncate to make keys more random
            // XXX: store key (lmax and maxi)
            uint32_t* r = lookup + lmax;
            unsigned j = *r;
            if (j == 1u) ++key_ct; // first instance of this key
            *r = b - s + 1;

            while ((c = gc(g)) != '>' && c != -1) {
                if (isspace(c)) {
                    fa->s_l = b - s;
                    assert((uint64_t)fa->s_l + fa_ent_max < (1ul << 32));
                    _buf_grow(fa->s, fa_ent_max);
                    b = s + fa->s_l;
                    continue;
                }
                if (i == width) {
                    pos += i;
                    i = 0;
                }
                c = b6(c);
                c = -isb6(c) & (c >> 1); // zero [min] for N - but also for A.

                if ((i & 3) == 0) *++b = '\0';
                *b |= c << ((i & 3) << 1);

                dna = (dna << 2) | c;
                rev = ((uint64_t)(c ^ 2) << KEYNT_TOP) | (rev >> 2);

                // select a strand dependent on the central Nt
                uint64_t key = dna & KEYNT_STRAND ? rev : dna;

                // with gray Nt code, we may select a more sequence specific section.
                uint64_t Ntgc = key ^ (key >> 2);

                if ((Ntgc > mx) || unlikely(Ntgc == mx && (dna & KEYNT_STRAND))) {
                    // a new maximum, keep it.
                    if (key == 0 && !(dna & KEYNT_STRAND))
                        break; // stretch of only 'N's or 'A's => skip

                    lmax = key;
                    maxi = i;
                    mx = Ntgc;

                } else if (i == maxi) {
                    // the max has left, we have to search for the next max
                    lmax = mx = Ntgc; // XXX: lmax converted below
                    j = maxi = i;

                    while ((j = (j - 1) % width) != i) {
                        Ntgc = rot[j] & ~KEYNT_STRAND;
                        Ntgc ^= Ntgc >> 2;
                        if (Ntgc > mx) {
                            mx = Ntgc;
                            lmax = rot[j];
                            maxi = j;
                        }
                    }
                } else {
                    rot[i] = key | (dna & KEYNT_STRAND);
                    ++i;
                    continue;
                }
                // XXX: store key (lmax and maxi)
                lmax = ((lmax >> 1) & ~HALF_KEYNT_MASK) | (lmax & HALF_KEYNT_MASK); // excise out 2nd cNt bit
                lmax &= KEYNT_TRUNC_MASK; // truncate to make keys more random
                r = lookup + lmax;

                j = *r;
                if (j == 1u) ++key_ct; // first element
                *r = fa->s_l; // exchange by this offset


                rot[i] = key | (dna & KEYNT_STRAND);
                ++i;
            }
        } while (c != '>' && c >= 0);
    } while (c >= 0);
    c = 0;
out: fa->key_ct = key_ct;
    return c;
}


