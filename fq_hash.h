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

#ifndef RK_FQ_HASH_H
#define RK_FQ_HASH_H

#include <stdint.h>
#include <assert.h>

#include <iostream>
#include <algorithm>

#include "khash.h"
#include "b6.h"

// maximum length of the read name
#define RK_FQ_MAX_NAME_ETC    (1ul << 14)

#define RK_FQ_KEYCOUNT_OFFSET 40

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

// exponentiate phred to allow summation
static unsigned phredtoe10[] = {
    1, 1, 1, 1, 2, 3, 3, 5, 6, 7, 10,
    12, 15, 19, 25, 31, 39, 50, 63, 79, 100,
    125, 158, 199, 251, 316, 398, 501, 630, 794, 1000,
    1258, 1584, 1995, 2511, 3162, 3981, 5011, 6309, 7943, 10000 };

using namespace std;

static inline unsigned bytes_for_bits(unsigned v)
{
    v >>= 3;
    v |= v >> 8;
    v |= v >> 16;
    return ++v;
}

struct fq_hash
{
    fq_hash(size_t rl, unsigned po) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + RK_FQ_MAX_NAME_ETC), nr(0), readlength(rl),
        phred_offset(po), keymax(0)
    {
        H = kh_init(UQCT);
        s = (uint8_t*)malloc(m);
        *s = '\0';
        max_bitm = bytes_for_bits(rl - 16);
        kroundup32(fq_ent_max);
    }

    /**
     * Store or update the hash with the given sequence. The key is the maximum of 16 Nts
     * surrounding a central-Nt, the strand decided by the cNts' 2nd bit. see encode_twisted().
     * The value is a 64 bit consisting of an offset in buffer `s' in
     * the low bits and a count of keys after RK_FQ_KEYCOUNT_OFFSET bits.
     *
     * At the offset, in the the hash value, buffer `s' holds four bytes - 1 bit, for the
     * offset of the max in the sequence,
     * and a `complement state' in the highest bit of the fourth byte. Another four bytes point
     * to the offset in buffer `s' with stats for the former read that had this same deviant,
     * or `1' if there were none.
     *
     * Finally the sequence and phred quality are encoded in the buffer, ended by 0xff, followed
     * by its read name - '\0' terminated.
     *
     * For each next read with the same deviant key, the hash values' offset in buffer `s' is
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

        // NOTE: the central Nt is not part of the key. A mismatch can occur at that Nt.
        khiter_t k = kh_get(UQCT, H, twisted);
        if (k == kh_end(H)) { /* new deviant */
            k = kh_put(UQCT, H, twisted, &is_err);
            if_ever ((is_err = (is_err < 0))) return;
            kh_val(H, k) = 1ul << RK_FQ_KEYCOUNT_OFFSET; // keycount
            i = 1; /* marks last element */

        } else {
//cerr << "Deviant 0x" << hex << deviant << dec << " already exists.\n";
            kh_val(H, k) += 1ul << RK_FQ_KEYCOUNT_OFFSET;
            i = kh_val(H, k) >> RK_FQ_KEYCOUNT_OFFSET;
            if (i > keymax) keymax = i;

            i = kh_val(H, k) & 0xffffffff; /* former value offset */
            kh_val(H, k) ^= i;
        }
        kh_val(H, k) |= l;

        for (unsigned j = 0; j != 4; ++j, i >>= 8)
            *++b = i & 0xff;

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
            khiter_t *ks, *kp;
            uint64_t *vp, *vs;
            uint64_t i = kh_size(H);
            uint8_t* dna = (uint8_t*)malloc((readlength<<1)+3);
            char seqrc[36];
            cout << hex;

            kp = ks = (khiter_t*)malloc(i * sizeof(khiter_t));

            for (khiter_t k = kh_begin(H); k != kh_end(H); ++k)
                if (kh_exist(H, k)) *kp++ = k;
            assert(kp != ks);
            cerr << "== " << i << " / " << nr << " sequences were unique," <<
                " the most frequent occured in " << keymax << " reads. ==\n";
            ++keymax;

            vs = (uint64_t*)malloc(keymax * sizeof(uint64_t));

            sort(ks, kp, *this); /* by keycount */

            do { /* loop over keys */
                i = kh_val(H, *--kp) & 0xffffffff;
                unsigned deviant = kh_key(H, *kp);
                decode_twisted(seqrc, s[i + max_bitm + 3] & 0x80, deviant);
                vp = vs;
                do { // loop over values belonging to key
                    *vp++ = i;
                    uint8_t* b = s + i + max_bitm + 7;
                    i = *b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b;
                } while (i != 1);

                sort(vs, vp, *this); // by offset/qualsum
                do {
                    uint8_t* b = s + *--vp + 8 + max_bitm;
                    unsigned len = decode_seqphred(b, dna);
                    cout << (char*)b + len + 1 << " 0x" << deviant << '\t' << seqrc << "\t";
                    b -= 8;
                    for (unsigned j = max_bitm; j != 0; --j) {
                        cout.width(2);
                        cout.fill('0');
                        cout << (unsigned)*--b;
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
        kh_destroy(UQCT, H);
    }
    inline bool operator() (khiter_t a, khiter_t b) /* needed by sort kp/ks: sort on keycounts*/
    {
        return kh_val(H, a) < kh_val(H, b);
    }

    inline bool operator() (uint64_t a, uint64_t b) /* needed by sort vs/vp: sort on orientation*/
    {
        uint8_t *sa = s + max_bitm + a + 3, *sb = s + max_bitm + b + 3;
        for (unsigned i = 0; i < 4; ++i, --sa, --sb)
            if (*sa != *sb) return *sa < *sb;
        return 0;
    }
    /*inline bool operator() (uint64_t a, uint64_t b)
    {
        uint8_t *sa = s + max_bitm + a, *sb = s + max_bitm + b;
        unsigned oa = 0u, ob = 0u;
        for (unsigned i = 0; i < 32; i += 8) {
            oa |= *sa++ << i;
            ob |= *sb++ << i;
        }
        oa &= ~(1u << 31);
        ob &= ~(1u << 31);
        return 0;
    }*/
    int is_err;

private:
    KHASH_INIT2(UQCT, kh_inline, uint32_t, uint64_t, 1, kh_int_hash_func, kh_int_hash_equal)

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
//cerr << "== " << sq << endl; //
        uint64_t lmax, rmax, c, d = 0;
        unsigned i = 0;
        do {
            if ((c = *sq++) == '\0') return 0;
            // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            c = b6(c);
            d = (d << 2) | (-isb6(c) & (c >> 1));

        } while (++i != 20);

        unsigned cNt = d & 0xc0000;     // excise out central Nt
        c = (0x40000 - 1) & d;          // mark set bits below it
        d ^= cNt ^ c ^ (c << 2);       // move lower Nts upward
        uint32_t t = cNt & 0x80000;  // 2nd bit of central Nt will determine the twist.
        uint64_t mCnt = cNt ^ t; //1st bit needs to be stored

        if ((c = *sq++) == '\0') return 0;
        c = b6(c);
        d |= -isb6(c) & (c >> 1);     // append next Nt
        *lmaxi = i;
        lmax = t ? d : (revseq(d) >> 24) ^ 0xaaaaaaaaaa; // truncated seq or revcmp according to 2nd bit of central Nt
        *m = 0;

        // if cNt bit is set: top part is seq, else revcmp
//cerr << "== Init: maxi " << *lmaxi << "\tmax: 0x" << hex << lmax << dec << endl;

        while((c = *sq++) != '\0') { // search for maximum from left
            c = b6(c);
            c = -isb6(c) & (c >> 1);

            d ^= cNt; // get next central Nt and put old one back;
            cNt ^= d & 0xc0000;
            t = cNt & 0x80000;
            d ^= cNt;

            d = (d << 2) | c ; // shift-insert next Nt.
            c = t ? (d & 0xffffffffff): (revseq(d) >> 24) ^ 0xaaaaaaaaaa;
            // need to truncate before maximize check.
            // rather than the max, we should search for the most distinct substring.
            if ((c > lmax) || ((c == lmax) && t)) { // yes, 9 x f : need to truncate before maximize check.
                lmax = c;
                *lmaxi = i;
//cerr << "== Mod: maxi " << *lmaxi << "\tmax: 0x" << hex << lmax << dec << endl;
                mCnt = cNt ^ t;
                *m |= 1u << (i & 7);
            }
            if ((++i & 7) == 0) *++m = 0u;
        }
        rmax = 0;
        while (--i != *lmaxi) { // search for maximum from right
            if ((i & 7) == 0) --m;
            c = *--sq;
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            d ^= cNt;
            cNt ^= d & 0xc0000;
            t = cNt & 0x80000;
            d ^= cNt;
            d = (c << 38) | (d >> 2);
            c = t ? d : (revseq(d) >> 16) ^ 0xaaaaaaaaaa;
            if ((c > rmax) || ((c == rmax) && !t)) {
                rmax = c;
                *m |= 1u << (i & 7);
            }
        }
//cerr << "== Res: maxi " << *lmaxi << "\tmax: 0x" << hex << lmax << dec << endl;
	*lmaxi |= -!mCnt & 0x80000000;
        return lmax & 0xffffffff; /* truncate on purpose */
    }
    void decode_twisted(char* seq, unsigned t, unsigned d)
    {
        char* revcmp = seq + 35;
        *revcmp = '\0';
	t = -!t & 2;

//cerr << "==t:" << t << endl;
        char c;

        for (unsigned i = 0; i != 16; ++i, ++seq, d >>= 2) {
            if (i == 10) { /* insert central Nt */
                *seq = b6(0 ^ t);
                *--revcmp = b6(4 ^ t);
                ++seq;
            }
            c = (d & 0x3) << 1;
            *seq = b6(4 ^ c);
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
            assert(*b <= 41);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' 0x30 is set and phred bits originally there shifted upwards
            // This preserves phred even for N, although this is probably unneccesary.
            *b |= isb6(c) ? (c << 5) : ((*b & 0x30) << 2) | 0x30;
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
            if ((c & 0x30) != 0x30) {
                *d = b6((c >> 5) & 0x6);
                *q = (c & 0x3f) + phred_offset;
            } else {
                *d = 'N';
                *q = (((c >> 2) & 0x30) | (c & 0xf)) + phred_offset;
            }
        }
        assert (len <= readlength);
        *d = *q = '\0';
        return len;
    }
    khash_t(UQCT) *H;
    size_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, keymax, max_bitm;
    uint8_t *s;
};


#endif //RK_FQ_HASH_H
