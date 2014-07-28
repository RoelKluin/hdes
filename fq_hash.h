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
#define KEY_MAX_CUTOFF 4
#define KEYNT_AC (1u << (KEY_WIDTH - 1))


#define CNtM        (KEYNT_AC - 1)
#define KEYCOUNT_OFFSET 31
#define KEY_BUFSZ (1u << (KEYCOUNT_OFFSET - 1))
#define BUFOFFSET_MASK ((1ul << KEYCOUNT_OFFSET) - 1)

//this must not exceed 32:
#define KEY_WIDTH (KEY_LENGTH + KEY_MAX_CUTOFF)

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

static inline void memset32(void* dest, uint32_t value, unsigned size)
{
  for(unsigned i = 0; i != size; i += 4)
    memcpy(((char*)dest) + i, &value, 4);
}
static inline void memset64(void* dest, uint64_t value, uint64_t size)
{
  for(uint64_t i = 0; i != size; i += 8)
    memcpy(((char*)dest) + i, &value, 8);
}

using namespace std;

static inline unsigned bytes_for_bits(unsigned v)
{
    v >>= 3;
    v |= v >> 8;
    v |= v >> 16;
    return ++v;
}
static inline bool key_cmp(uint64_t a, uint64_t b) /* needed by sort kp/ks: sort on keycounts*/
{
    return (a & ~BUFOFFSET_MASK) < (b & ~BUFOFFSET_MASK);
    //return (a & BUFOFFSET_MASK) < (b & BUFOFFSET_MASK);
}

struct fq_arr
{
    fq_arr(size_t rl, unsigned po) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + RK_FQ_MAX_NAME_ETC), nr(0), readlength(rl),
        phred_offset(po)
    {
        s = (uint8_t*)malloc(m);
        *s = '\0';
        rl -= KEY_WIDTH;
        key_ct = 0u;
        max_bitm = bytes_for_bits(rl);
        unsigned looklen = KEY_BUFSZ;
        look = (uint64_t*)malloc(sizeof(uint64_t)*looklen);
        assert(look != NULL);
        memset64(look, 1ul, sizeof(uint64_t)*looklen);
        kroundup32(fq_ent_max);
    }

    /**
     * Store or update the entry with the given sequence. The key is the maximum of 16 Nts
     * surrounding a central-Nt, the strand decided by the cNts' 2nd bit. see encode_twisted().
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
        twisted = encode_twisted(seq, b, &i); // remove low bit, separate deviant
        if_ever (i == 0) return;           // XXX: means sequence too short. skipped.
        b += max_bitm - 1;

        for (unsigned j = 0; j != 4; ++j, i >>= 8) // insert offset
            *++b = i & 0xff;

        // TODO: do without KEYCOUNT (count per key - keep key_ct, i.e. tot key cpunt)
        size_t* r = &look[twisted];
        size_t t = *r & BUFOFFSET_MASK;

        if (t == 1ul) key_ct++; // first element
        assert (l < KEY_BUFSZ);
        *r ^= t ^ l; // exchange by this offset
        *r += 1ul << KEYCOUNT_OFFSET;

        for (unsigned j = 0; j != 4; ++j, t >>= 8)
            *++b = t & 0xff; // move former offset, or last element mark (1u) here

        b = encode_seqphred(qual, seq, b);
        while ((*++b = *name) != '\0') ++name;
        l = ++b - s;
    }

    /**
     * print the fastq entries in a sorted order. Fastq entries for deviant keys
     * that occur the most first. Sorted on the phred sum for identical counts.
     * for each deviant the reads are then sorted on the offset of the deviant
     * key in the read. This occurs for both strands.
     */
    void dump() { /* must always be called to clean up */
        if (not is_err) {
            unsigned *ks, *kp;
            uint64_t *vp, *vs;
            uint64_t i = 0u, keymax = 1ul << KEYCOUNT_OFFSET;
            uint8_t* dna = (uint8_t*)malloc((readlength<<1)+3);
            char seqrc[(KEY_LENGTH+1)<<1];
            cout << hex;

            kp = ks = (unsigned*)malloc(key_ct * sizeof(unsigned));

            unsigned uniq = 0;
            uint64_t j = KEY_BUFSZ - 1;
            do {
                if (look[j] != 1ul) { // has keys
                    *kp++ = j;
                    i = look[j];
                    i &= ~BUFOFFSET_MASK;
                    if (i == (1ul << KEYCOUNT_OFFSET)) ++uniq;
                    else if (i > keymax) keymax = i;
                }
            } while (--j != 0);
            keymax >>= KEYCOUNT_OFFSET;

            assert(kp != ks);
            cerr << "== " << key_ct << " / " << nr << " keys per sequences," <<
                " of which " << uniq << " were unique." <<
                " the most frequent occured in " << keymax << " reads. ==\n";
            ++keymax;

            // FIXME: why is keymax not sufficient (valgrind on very large fq)? 
            // maybe due to key collision?
            // CTCTGTTGTGTCTGATT|AATCAGACACAACAGAG
            // CTCTGTGGTGTCTGATT|AATCAGACACCACAGAG
            //       ^                     ^
            vs = (uint64_t*)malloc((keymax << 1) * sizeof(uint64_t));

            sort(ks, kp, key_cmp); /* by keycount */
            // TODO merge two loops, already [sorted on key]

            do { /* loop over keys */
                j = *--kp;
                decode_twisted(seqrc, j);
                i = look[j] & BUFOFFSET_MASK;
                assert(i != 1ul);
                vp = vs;
                do { // loop over values belonging to key
                    *vp++ = i; // XXX: valgrind: invalid write here?
                    uint8_t* b = s + i + max_bitm + 7;
                    i = *b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b;
                    assert (i < KEY_BUFSZ);
                } while (i != 1ul);

                sort(vs, vp, *this); // by offset/qualsum
                do {
                    uint8_t* b = s + *--vp + 8 + max_bitm;
                    unsigned len = decode_seqphred(b, dna);
                    cout << (char*)b + len + 1 << " " <<
                        " 0x" << j << '\t' << seqrc << "\t";
                    b -= 8 + max_bitm;
                    for (unsigned j = max_bitm; j != 0; --j) {
                        uint8_t c = BitReverseTable256[*b++];
                        cout << (unsigned)((c & 0xf0) >> 4);
                        cout << (unsigned)(c & 0xf);
                    }
                    cout << endl << dna << "\n+\n";
                    cout << dna + readlength + 1 << "\n";
                } while (vp != vs);
            } while (kp != ks);
            free(dna);
            free(vs);
            free(ks);
        }
	free(s);
        free(look);
    }

    inline bool operator() (uint64_t a, uint64_t b) /* needed by sort vs/vp: sort on orientation*/
    {
        uint8_t *sa = s + max_bitm + a + 3, *sb = s + max_bitm + b + 3;
        for (unsigned i = 0; i < 4; ++i, --sa, --sb)
            if (*sa != *sb) return *sa < *sb;
        return 0;
    }
    int is_err;

private:

    /**
     * The sequence is converted to revcmp according to the central Nt, which is not part of the
     * sequence.
     *
     * The deviant part is insensitive to the strand. The sequence part is converted to the other
     * strand according to the 2nd bit of the central Nt. This ensures also the sequence part in
     * the seqdeviant is the same, regardless of the original strand.
     *
     * For all promising 17 Nt windows fitting our sequences' readlength, the seqdev is calculated.
     * The maximum thereof is kept as well as its offset and the first bit of the central Nt.
     *
     */
    uint32_t encode_twisted(const uint8_t* sq, uint8_t* m, unsigned* lmaxi)
    {
        uint64_t lmax, maxv, c, d = 0, r = 0;
        int i = -KEY_WIDTH;
        do {
            if ((c = *sq++) == '\0') return 0;
            // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            d = (d << 2) | c;
            r = (c << KEY_TOP_NT) | (r >> 2);

        } while (++i < 0);

        r ^= KEY_WIDTH_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement
        // seq or revcmp according to 2nd bit of central Nt
        c = (d & KEYNT_STRAND) ? d : r;
        lmax = c;
        maxv = c ^ (c >> 2); // maximize on gray Nt sequence
        *lmaxi = i;
        *m = 1;

        // if cNt bit is set: top part is seq, else revcmp

        while((c = *sq++) != '\0') { // search for maximum from left
            if ((++i & 7) == 0) *++m = 0u;
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            d = ((d << 2) & KEY_WIDTH_MASK & 0xffffffffffffffff) | c;
            r = ((c ^ 2) << KEY_TOP_NT) | (r >> 2);
            c = (d & KEYNT_STRAND) ? d : r;
            uint64_t Ntgc = c ^ (c >> 2);

            // need to truncate before maximize check.
            // rather than the max, we should search for the most distinct substring.
            if ((Ntgc > maxv) || unlikely((Ntgc == maxv) && (d & KEYNT_STRAND))) {
                lmax = c;
                maxv = Ntgc;
                *lmaxi = i;
                *m |= 1u << (i & 7);
            }
        }
        *m |= 1u << (i & 7);
        sq -= KEY_WIDTH + 2;
        c = (d & KEYNT_STRAND) ? d : r;
        maxv = c ^ (c >> 2);

        while ((unsigned)--i > *lmaxi) { // search for maximum from right
            if ((i & 7) == 0) --m;
            c = *sq--;
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            d = (c << KEY_TOP_NT) | (d >> 2); // Note: walking back
            r = ((r << 2) & KEY_WIDTH_MASK & 0xffffffffffffffff) | (c ^ 2);
            c = (d & KEYNT_STRAND) ? d : r;
            uint64_t Ntgc = c ^ (c >> 2);

            if ((Ntgc > maxv) || unlikely((Ntgc == maxv) && !(d & KEYNT_STRAND))) {
                maxv = Ntgc;
                *m |= 1u << (i & 7);
            }
        }
        lmax = ((lmax >> 1) & ~HALF_KEY_WIDTH_MASK) | (lmax & HALF_KEY_WIDTH_MASK);     // excise out central Nt
        return lmax & (KEY_BUFSZ - 1u); /* truncate on purpose */
    }
    void decode_twisted(char* seq, unsigned d)
    {
        char* revcmp = seq + (KEY_LENGTH << 1) + 1;
        *revcmp = '\0';

        char c;

        for (unsigned i = 0; i != KEY_LENGTH - 1; ++i, ++seq, d >>= 2) {
            if (i == ((KEY_WIDTH - 1) >> 1)) { /* insert central Nt */
                c = (d & 1) << 1;
                *seq = b6(c);
                *--revcmp = b6(c ^ 4);
                ++seq;
                d >>= 1;
            }
            c = (d & 0x3) << 1;
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
        while ((*++b = *qual++) != '\0') {
            *b -= phred_offset;
            assert(*b <= 50);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' phred should always be less than 3
            // max phred for base (G) is 60 using this scheme
            if (isb6(c)) *b = (*b | (c << 5)) + 3;
        }
        *b = 0xff;
        return b;
    }
    /**
     * decode the sequence and phreds.
     */
    unsigned decode_seqphred(uint8_t* b, uint8_t* d)
    {
        unsigned len;
        uint8_t *q = d + readlength + 1;

        for(len = 0; *b != 0xff; ++len, ++b, ++d, ++q) {
            unsigned c = *b;
            if (c >= 3) {
                c -= 3;
                *d = b6((c >> 5) & 0x6);
                *q = (c & 0x3f) + phred_offset;
            } else {
                *d = 'N';
                *q = c + phred_offset;
            }
        }
        assert (len <= readlength);
        *d = *q = '\0';
        return len;
    }
    size_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, max_bitm, key_ct;
    uint8_t *s;
    uint64_t* look;
};


#endif //FQ_ARR_H

