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

#define MAX_NAME_ETC    4096

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

static unsigned phredtoe10[] = {
    1, 1, 1, 1, 2, 3, 3, 5, 6, 7, 10,
    12, 15, 19, 25, 31, 39, 50, 63, 79, 100,
    125, 158, 199, 251, 316, 398, 501, 630, 794, 1000,
    1258, 1584, 1995, 2511, 3162, 3981, 5011, 6309, 7943, 10000 };

using namespace std;

static inline uint32_t revseq(uint32_t x)
{
    uint32_t m = 0x33333333;
    x = ((x & m) << 2) | ((x & ~m) >> 2);
    m = 0x0f0f0f0f;
    x = ((x & m) << 4) | ((x & ~m) >> 4);
    m = 0x00ff00ff;
    x = ((x & m) << 8) | ((x & ~m) >> 8);
    m = 0x0000ffff;
    return ((x & m) << 16) | ((x & (m << 16)) >> 16);
}

struct fq_hash
{
    fq_hash(size_t rl, unsigned po) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + MAX_NAME_ETC), nr(0), readlength(rl), phred_offset(po), keymax(0)
    {
        H = kh_init(UQCT);
        s = (uint8_t*)malloc(m);
        dna = (uint8_t*)malloc((rl<<1)+2);
        *s = '\0';
        kroundup32(fq_ent_max);
    }
    /*
     * KC-SL-key_1 PS-KO-val_1 KC-SL-key_2] PS-KO-val_2 ... KC-SL-key_n-1 PS-KO-val_n-1 PS-*VO_1*-val_n]
     *
     * KC:keycount, SL: seqlength. PS:phred sum, KO key offset, VO_1: Value 1 offset.
     *
     * val: phred + name, key: dna in 2bit. valn has existing key1 (same for val1)
     */
    void put(const uint8_t* seq, const uint8_t* qual, const uint8_t* name)
    {
        if_ever (is_err) return;
        ++nr;
        if (l + fq_ent_max >= m) { // grow buffer if insufficient space for another read
            m <<= 1;
            s = (uint8_t*)realloc(s, m);
        }
        unsigned i = 0;
        // get twisted: value that is the same for seq and revcmp
        unsigned twisted = encode_twisted(seq, &i); // remove low bit, separate deviant
        if_ever (i == 0) return; //sequence too short.

        uint8_t *b = s + l + 4; // four bytes for phred sum
        *b = (i & 0xff); i >>= 8; // insert offset / twisted state
        *++b = (i & 0xff); i >>= 8;
        *++b = (i & 0xff); i >>= 8;
        *++b = (i & 0xff);

        // NOTE: the central Nt is not part of the key. A mismatch can occur at that Nt.
        khiter_t k = kh_get(UQCT, H, twisted);
        if (k == kh_end(H)) { /* new deviant */
            k = kh_put(UQCT, H, twisted, &is_err);
            if_ever ((is_err = (is_err < 0))) return;
            //TODO: could increase pot. buffersize at the expense of max keycounts
            kh_val(H, k) = 1ul << 32; // keycount
            i = 1; /* marks last element */

        } else {
//cerr << "Deviant 0x" << hex << deviant << dec << " already exists.\n";
            kh_val(H, k) += 1ul << 32;
            i = kh_val(H, k) >> 32;
            if (i > keymax) keymax = i;

            i = kh_val(H, k) & 0xffffffff; /* former value offset */
            kh_val(H, k) ^= i;
        }
        // TODO: offset in seq and devbit
        kh_val(H, k) |= l;

        /* former value offset */
        for (unsigned j = 0; j != 4; ++j) {
            *++b = i & 0xff;
            i >>= 8;
        }
        b = encode_seqphred(qual, seq, b);
        while ((*++b = *name) != '\0') ++name;
        l = ++b - s;
    }

    void dump() { /* must always be called to clean up */
        if (not is_err) {
            khiter_t *ks, *kp;
            uint64_t *vp, *vs;
            uint64_t i = kh_size(H);
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
                decode_twisted(seqrc, s[i + 4] & 1, deviant);
                vp = vs;
                do { // loop over values belonging to key
                    *vp++ = i;
                    uint8_t* b = s + i + 11;
                    i = *b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b;
                } while (i != 1);

                sort(vs, vp, *this); // by offset/qualsum
                do {
                    unsigned len = decode_seqphred(s + *--vp + 12);
                    cout << (char*)s + *vp + 13 + len << " 0x" << deviant << '\t' << seqrc << endl;
                    cout << dna << "\n+\n";
                    cout << dna + readlength + 1 << "\n";
                } while (vp != vs);
            } while (kp != ks);
            free(vs);
            free(ks);
        }
        free(dna);
	free(s);
        kh_destroy(UQCT, H);
    }
    inline bool operator() (khiter_t a, khiter_t b) /* needed by sort kp/ks: sort on keycounts*/
    {
        return kh_val(H, a) < kh_val(H, b);
    }

    inline bool operator() (uint64_t a, uint64_t b) /* needed by sort vs/vp: sort on orientation/phred sum*/
    {
        uint8_t *sa = s + a + 7, *sb = s + b + 7;
        for (unsigned i = 0; i < 8; ++i, --sa, --sb)
            if (*sa != *sb) return *sa < *sb;
        return 0;
    }
    int is_err;

private:
    KHASH_INIT2(UQCT, kh_inline, uint32_t, uint64_t, 1, kh_int_hash_func, kh_int_hash_equal)
    unsigned get_twisted(unsigned rev, unsigned sq)
    {
	rev = rev ^ 0xaaaaaaaau;
	unsigned dev = sq ^ rev;

	unsigned m = dev & -dev; // leftmost deviant bit. NB: palindromes have none
	sq &= m;                 // test whether it's set for first sequence
	rev ^= (sq ^ -sq) & dev; // flip revcmp deviant above devbit, if set

	m -= !!m; // mask below devbit, if any - none for palindromes
	m &= rev; // revnxt below devbit, if any
	rev ^= m ^ (m << 1) ^ sq ^ !sq; // swap devbit to lowest and flip
	return ((rev & 0xffff) << 16) | (dev & 0xffff);
    }
    /* returns no. characters processed. */
    unsigned encode_twisted(const uint8_t* sq, unsigned* maxi)
    {
//cerr << "== " << sq << endl; //
        unsigned i = 0;
        uint32_t d = 0;
        uint32_t c;
        do {
            if ((c = *sq++) == '\0') return 0;
            // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            c = b6(c);
            d = (d << 2) | (-isb6(c) & (c >> 1));

        } while (++i != 16);

        unsigned cNt = d & 0xc000;     // excise out central Nt
        c = (0x4000 - 1) & d;          // mark set bits below it
        uint32_t t = !(cNt & 0x8000);  // first bit of central Nt will determine the orientation.

        d ^= cNt ^ c ^ (c << 2);       // move lower Nts upward

        if ((c = *sq++) == '\0') return 0;
        c = b6(c);
        d |= -isb6(c) & (c >> 1);     // append next Nt
        *maxi = (i << 1) | ((cNt & 0x4000) == 0); // init: offset and cNt bit, orientation.
        t = -t | 0xffff;

        t &= d ^ revseq(d) ^ 0xaaaaaaaa; // create deviant half or whole - dependent on sinbit
        unsigned max = t ^ (d & 0xffff0000); // high bits become seq or revcmp accordingly

        // if cNt bit is set: top part is seq, else revcmp
//cerr << "== Init: maxi " << *maxi << "\tmax: 0x" << hex << max << dec << endl;

        do { // search for maximum deviant
//cerr << "==:" << i << endl;
            if ((c = *sq++) == '\0') break;
            c = b6(c);
            c = -isb6(c) & (c >> 1);

            d ^= cNt; // get next central Nt and put old one back;
            cNt ^= d & 0xc000;
            t = !(cNt & 0x8000);
            d ^= cNt;

            d = ((d & 0x3fffffff) << 2) | c ; // shift-insert next Nt.

            if ((t == 0) && ((d | 0xffff) < max)) continue;

            c = (-t | 0xffff) & (d ^ revseq(d) ^ 0xaaaaaaaa);
            c ^= d & 0xffff0000;

            // rather than the max, we should search for the most distinct substring.
            if (c > max) {
                max = c;
                *maxi = (i << 1) | ((cNt & 0x4000) == 0);
//cerr << "== Mod: maxi " << *maxi << "\tmax: 0x" << hex << max << dec << endl;
            }
        } while (++i != readlength);
//cerr << "== Res: maxi " << *maxi << "\tmax: 0x" << hex << max << dec << endl;
        return max;
    }
    void decode_twisted(char* seq, unsigned t, unsigned d)
    {
        char* revcmp = seq + 35;
        *revcmp = '\0';

        unsigned i;
//cerr << "==t:" << t << endl;

        d ^= (revseq(d) & 0xffff) ^ 0xaaaa;
        char c;

        for (i = 0; i != 8; ++i, ++seq, d >>= 2) {
            c = (d & 0x3) << 1;
            *seq = b6(4 ^ c);
            *--revcmp = b6(c);
        }
        t <<= 1;
        *seq = b6(2 ^ t);
        *--revcmp = b6(6 ^ t);

        for (++seq; i != 16; ++i, ++seq, d >>= 2) {
            c = (d & 0x3) << 1;
            *seq = b6(4 ^ c);
            *--revcmp = b6(c);;
        }
        *seq = '|';
        assert(revcmp == seq + 1);
    }

    uint8_t* encode_seqphred(const uint8_t* qual, const uint8_t* seq, uint8_t* b)
    {
        unsigned i = 0;
        uint8_t* p = b - 11;
        while ((*++b = *qual++) != '\0') {
            *b -= phred_offset;
            assert(*b <= 40);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' 0x30 is set and phred bits originally there shifted upwards
            // This preserves phred even for N, although this is probably unneccesary.
            *b |= isb6(c) ? (c << 5) : ((*b & 0x30) << 2) | 0x30;
        }
        *b = 0xff;
        *p = (i & 0xff); i >>= 8; // insert phred sum
        *++p = (i & 0xff); i >>= 8;
        *++p = (i & 0xff); i >>= 8;
        *++p = (i & 0xff);
        return b;
    }
    unsigned decode_seqphred(uint8_t* b)
    {
        unsigned len;
        uint8_t *d = dna, *q = dna + readlength + 1;

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
        *d = *q = '\0';
        return len;
    }
    khash_t(UQCT) *H;
    size_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, keymax;
    uint8_t *s;
    uint8_t* dna;
};


#endif //RK_FQ_HASH_H

