/*
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel KLuin, 
 */

#ifndef FQ_ARR_H
#define FQ_ARR_H

#include <stdint.h>
#include <assert.h>

#include <iostream>
#include <algorithm>
#include <string.h>
#include "b6.h"
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


// maximum length of the read name
#define RK_FQ_MAX_NAME_ETC    (1ul << 14)

// KEY_LENGTH must be odd - or 2nd bit of central Nt is not always flipped in its
// complement - the alternative is the twisted halfdev conversion, but this is cheaper
#define KEY_LENGTH 15
#define KEY_CUTOFF 4
//this must not exceed 32:
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF)
#define KEYNT_AC (1u << (KEY_WIDTH - 1))


#define CNtM        (KEYNT_AC - 1)
#define KEYCOUNT_OFFSET 31
#define KEY_BUFSZ (1u << (KEYCOUNT_OFFSET - 1))
#define BUFOFFSET_MASK ((1ul << KEYCOUNT_OFFSET) - 1)

#define KEYNT_STRAND (1u << KEY_WIDTH)
#define KEY_TOP_NT ((KEY_WIDTH - 1) * 2)

#define KEY_WIDTH_MASK ((1ul << KEY_WIDTH * 2) - 1ul)
#define HALF_KEY_WIDTH_MASK (KEYNT_STRAND - 1u)

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

// exponentiate phred to allow summation
// as.integer(round(10^((0:63)/10)))
static unsigned phredtoe10[] = {
    1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10, 13, 16, 20,
    25, 32, 40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501,
    631, 794, 1000, 1259, 1585, 1995, 2512, 3162, 3981, 5012, 6310, 7943, 10000, 12589,
    15849, 19953, 25119, 31623, 39811, 50119, 63096, 79433, 100000, 125893, 158489, 199526,
    251189, 316228, 398107, 501187, 630957, 794328, 1000000, 1258925, 1584893, 1995262
};

template<typename T>
static inline void memset_T(void* dest, T value, T size)
{
  for(T i = 0; i != size; i += sizeof(T))
    memcpy(((char*)dest) + i, &value, sizeof(T));
}

using namespace std;

struct fq_arr
{
    fq_arr(size_t rl, unsigned po) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + RK_FQ_MAX_NAME_ETC), nr(0), readlength(rl),
        phred_offset(po)
    {
        s = (uint8_t*)malloc(m);
        *s = '\0';
        key_ct = 0u;
        unsigned looklen = KEY_BUFSZ;
        look = (uint64_t*)malloc(sizeof(uint64_t)*looklen);
        assert(look != NULL);
        memset_T<uint64_t>(look, 1ul, sizeof(uint64_t)*looklen);
        kroundup32(fq_ent_max);
    }

    /**
     * Store or update the entry with the given sequence. The key is the maximum of 16 Nts
     * surrounding a central-Nt, the strand decided by the cNts' 2nd bit. see maximize_key().
     * The value is a 64 bit consisting of an offset in buffer `s' in
     * the low bits and a count of keys after KEYCOUNT_OFFSET bits.
     *
     * At the offset, buffer `s' holds four bytes - 1 bit, for the
     * offset of the max in the sequence,
     * and a `complement state' in the highest bit of the fourth byte. Another four bytes point
     * to the offset in buffer `s' with stats for the former read that had this same deviant,
     * or `1' if there were none.
     *
     * Finally the sequence and phred quality are encoded in the buffer, ended by 0xff, followed
     * by its read name - '\0' terminated.
     *
     * For each next read with the same deviant key, the values' offset in buffer `s' is
     * updated to point to its last entry, until `1' which instead shows it was the last entry.
     */
    void put(const uint8_t* seq, const uint8_t* qual, const uint8_t* name)
    {
        if_ever (is_err) return;
        assert (l != 1);
        ++nr;
        if (l + fq_ent_max >= m) { // grow buffer if insufficient space for another read
            m <<= 1;
            s = (uint8_t*)realloc(s, m);
        }
        unsigned twisted, i = 0;
        uint8_t *b = s + l; // bytes to store max locations
        // get twisted: value that is the same for seq and revcmp
        i = init_key(seq);
        if_ever (i == 0) return;           // XXX: means sequence too short. skipped.
        twisted = maximize_key(seq + KEY_WIDTH, &i);

        for (unsigned j = 0; j != 4; ++j, i >>= 8) // insert key offset
            *b++ = i & 0xff;

        // TODO: do without KEYCOUNT (count per key - keep key_ct, i.e. tot key cpunt)
        size_t* r = &look[twisted];
        size_t t = *r & BUFOFFSET_MASK;

        if (t == 1ul) key_ct++; // first element
        assert (l < KEY_BUFSZ);
        *r ^= t ^ l; // exchange by this offset
        *r += 1ul << KEYCOUNT_OFFSET;

        for (unsigned j = 0; j != 4; ++j, t >>= 8)
            *b++ = t & 0xff; // move former buffer offset, or last element mark (1u) here

        b = encode_seqphred(qual, seq, b);
        while ((*b++ = *name) != '\0') ++name;
        l = b - s;
    }

    /**
     * print the fastq entries in a sorted order. Fastq entries for deviant keys
     * that occur the most first. Sorted on the phred sum for identical counts.
     * for each deviant the reads are then sorted on the offset of the deviant
     * key in the read. This occurs for both strands.
     */
    void dump() { /* must always be called to clean up */
        if (not is_err) {
            uint64_t i = 0u;
            uint8_t* d = (uint8_t*)malloc((readlength<<1)+3);
            char seqrc[(KEY_LENGTH+1)<<1];
            cout << hex;

            uint64_t j = KEY_BUFSZ - 1;
            do {
                i = look[j] & BUFOFFSET_MASK;
                if (i != 1ul) {
                    decode_key(seqrc, j);
                    do { // has keys
                        uint8_t* b = s + i + 7;
                        i = *b; i <<= 8;
                        i |= *--b; i <<= 8;
                        i |= *--b; i <<= 8;
                        i |= *--b;
                        assert (i < KEY_BUFSZ);
                        b += 4;
                        unsigned len = decode_seqphred(b, d);
                        cout << (char*)b + len + 1 << " " <<
                            " 0x" << j << '\t' << seqrc << endl;

                        cout << d << "\n+\n" << d + readlength + 1 << "\n";
                    } while (i != 1ul);
                }
            } while (--j != 0);
            free(d);
        }
	free(s);
        free(look);
    }
    int is_err;

private:

    uint32_t init_key(const uint8_t* sq)
    {
        uint64_t c;
        dna = rev = 0ul;

        for (rot = KEY_WIDTH; rot != 0; --rot) {
            if ((c = *sq++) == '\0') return 0;
            // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            dna = (dna << 2) | c;
            rev = (c << KEY_TOP_NT) | (rev >> 2);
        }
        rev ^= KEY_WIDTH_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement
        return KEY_WIDTH;
    }
    /**
     * The key is the sequence or revcmp dependent on the 2nd bit of the central Nt. As every
     * other 2nd bit of each twobit, it is inverted for the reverse complement. Since this is
     * the central Nt - odd KEY_LENGTH - the position is the same for the reverse complement.
     *
     * For all windows of KEY_LENGTH plus KEY_CUTOFF, the key value is calculated. The maximum
     * thereof is kept as well as its offset.
     *
     */
    uint32_t maximize_key(const uint8_t* sq, unsigned* lmaxi)
    {
        uint64_t lmax, maxv, c;

        // seq or revcmp according to 2nd bit of central Nt
        c = (dna & KEYNT_STRAND) ? dna : rev;
        lmax = c;
        maxv = c ^ (c >> 2); // maximize on gray Nt sequence
        *lmaxi = rot;

        // if cNt bit is set: top part is seq, else revcmp

        while((c = *sq++) != '\0') { // search for maximum from left
            ++rot;
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            dna = ((dna << 2) & KEY_WIDTH_MASK & 0xffffffffffffffff) | c;
            rev = ((c ^ 2) << KEY_TOP_NT) | (rev >> 2);
            c = (dna & KEYNT_STRAND) ? dna : rev;
            uint64_t Ntgc = c ^ (c >> 2);

            // need to truncate before maximize check.
            // rather than the max, we should search for the most distinct substring.
            if ((Ntgc > maxv) || unlikely((Ntgc == maxv) && (dna & KEYNT_STRAND))) {
                lmax = c;
                maxv = Ntgc;
                *lmaxi = rot;
            }
        }
        *lmaxi += KEY_WIDTH;//
        lmax = ((lmax >> 1) & ~HALF_KEY_WIDTH_MASK) | (lmax & HALF_KEY_WIDTH_MASK);     // excise out central Nt
        return lmax & (KEY_BUFSZ - 1u); /* truncate on purpose */
    }
    void decode_key(char* seq, unsigned dna)
    {
        char* revcmp = seq + (KEY_LENGTH << 1) + 1;
        *revcmp = '\0';

        char c;

        for (unsigned i = 0; i != KEY_LENGTH - 1; ++i, ++seq, dna >>= 2) {
            if (i == ((KEY_WIDTH - 1) >> 1)) { /* insert central Nt */
                c = (dna & 1) << 1;
                *seq = b6(c);
                *--revcmp = b6(c ^ 4);
                ++seq;
                dna >>= 1;
            }
            c = (dna & 0x3) << 1;
            *seq = b6(c ^ 4);
            *--revcmp = b6(c);
        }
        *seq = '|';
        assert(revcmp == seq + 1);
    }

    /**
     * put sequence and phred qualities in one char foreach Nt. This includes the deviant
     * sequence section currently, but this may change in the future - it is not needed.
     */
    uint8_t* encode_seqphred(const uint8_t* qual, const uint8_t* seq, uint8_t* b)
    {
        unsigned i = 0;
        while ((*b = *qual++) != '\0') {
            *b -= phred_offset;
            assert(*b <= 50);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' phred should always be less than 3
            // max phred for base (G) is 60 using this scheme
            if (isb6(c)) *b = (*b | (c << 5)) + 3;
            ++b;
        }
        *b++ = 0xff;
        return b;
    }
    /**
     * decode the sequence and phreds.
     */
    unsigned decode_seqphred(uint8_t* b, uint8_t* d)
    {
        unsigned len;
        uint8_t *r = d + readlength + 1;

        for(len = 0; *b != 0xff; ++len, ++b, ++d, ++r) {
            unsigned c = *b;
            if (c >= 3) {
                c -= 3;
                *d = b6((c >> 5) & 0x6);
                *r = (c & 0x3f) + phred_offset;
            } else {
                *d = 'N';
                *r = c + phred_offset;
            }
        }
        assert (len <= readlength);
        *d = *r = '\0';
        return len;
    }
    size_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, key_ct, rot;
    uint8_t *s;
    uint64_t* look, dna, rev;
};


#endif //FQ_ARR_H

