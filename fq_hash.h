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

#define MAX_NAME_ETC    1024

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

struct fq_hash
{
    fq_hash(size_t rl) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + MAX_NAME_ETC), nr(0), readlength(rl), phred_offset(33), keymax(0)
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
        if (l + fq_ent_max >= m) {
            m <<= 1;
            s = (uint8_t*)realloc(s, m);
        }
        uint8_t *p, *b = s + l;
        size_t i = encode_2bit(seq);

        khiter_t k = str_kh_get(H, (const char*)dna);
        if (k == kh_end(H)) { /* new key */
            for (p = b + 8; b != p; ++b) *b = '\0'; /* keycount(5)/length(3) */
            unsigned c = i;
            *--p = (i & 0xff); i >>= 8; /* seqlength */
            *--p = (i & 0xff); i >>= 8;
            *--p = (i & 0xff);
            i = p - s - 5; /* key offset, (denotes last value) */
            p = dna;
            while (c--) *b++ = *p++;
            ++b;
            k = kh_put(UQCT, H, i, &is_err);
            if_ever ((is_err = (is_err < 0))) return;

        } else {
            p = s + kh_key(H, k) + 3;

            /* increment keycount */
            i = ++(*p);
            i |= (*--p += !i) << 8;
            i |= (*--p += !i) << 16;
            i |= (*--p += !i) << 24;
            if (i > keymax) keymax = i;
            i = kh_val(H, k); /* former value offset */
        }
        // 4 bytes reserved for phred sum
        p = b;
        b += 3;
        /* put former value offset here */
        for (unsigned j = 0; j != 4; ++j) {
            *++b = i & 0xff;
            i >>= 8;
        }
        b = encode_phredqual(qual, seq, b);
        while ((*++b = *name) != '\0') ++name;
        l = ++b - s;
        kh_val(H, k) = p - s;
    }

    void dump() { /* must always be called to clean up */
        if (not is_err) {
            khiter_t *ks, *kp;
            size_t *vp, *vs;
            s[l] = '\0'; // end last value
            size_t i = kh_size(H);

            kp = ks = (khiter_t*)malloc(kh_size(H) * sizeof(khiter_t));
            for (khiter_t k = kh_begin(H); k != kh_end(H); ++k)
                if (kh_exist(H, k)) *kp++ = k;
            ++keymax;
            cerr << "== " << i << " / " << nr << " sequences were unique," <<
                " the most frequent occured in " << keymax << " reads. ==\n";

            vs = (size_t*)malloc(keymax * sizeof(size_t));

            sort(ks, kp, *this); /* by keycount */

            while (kp != ks) { /* loop over keys */
                size_t sq = kh_key(H, *--kp);
                i = kh_val(H, *kp);
                uint8_t* p = s + sq + 5;
                unsigned len = *p; len <<= 8;
                len |= *++p; len <<= 8;
                len |= *++p;
                assert(len < (1 << 10)); //


                vp = vs;
                do { // loop over values belonging to key
                    *vp++ = i;
                    uint8_t* b = s + i + 8;
                    i = *--b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b; i <<= 8;
                    i |= *--b;

                } while (i != sq);

                sort(vs, vp, *this); // by qualsum
                while (vp-- != vs) {
                    decode_phredqual(s + *vp + 8, len);
                    cout << (char*)s + *vp + 9 + len << "\n";
                    cout << dna << "\n+\n";
                    cout << dna + len + 1 << "\n";
                }
            }
            free(vs);
            free(ks);
        }
        free(dna);
	free(s);
        kh_destroy(UQCT, H);
    }
    inline bool operator() (khiter_t a, khiter_t b) /* needed by sort kp/ks.*/
    {
        uint8_t *sa = s + kh_key(H, a);
        uint8_t *sb = s + kh_key(H, b);
        while (*sa == *sb) ++sa, ++sb; /* No need to check length, no 2 keys are the same.*/
        return *++sa < *++sb;
    }

    inline bool operator() (size_t a, size_t b) /* needed by sort vs/vp*/
    {
        uint8_t *sa = s + a, *sb = s + b;
        if (*sa != *sb) return *sa < *sb;
        if (*++sa != *++sb) return *sa < *sb;
        if (*++sa != *++sb) return *sa < *sb;
        return *++sa < *++sb;
    }
    int is_err;

private:
#define my_str_hash_equal(a, b) (strcmp((char*)s + a + 8, (char*)s + b + 8) == 0)
#define my_str_hash_func(key) __ac_X31_hash_string((char*)s + key + 8)
    KHASH_INIT2(UQCT, static kh_inline, uint32_t, uint32_t, 1, my_str_hash_func, my_str_hash_equal)
    static kh_inline khint_t str_kh_get(const kh_UQCT_t *h, const char* key)
    {
        if (h->n_buckets) {
            khint_t inc, k, i, last, mask;
            mask = h->n_buckets - 1;
            k = __ac_X31_hash_string(key); i = k & mask;
            inc = __ac_inc(k, mask); last = i; // inc==1 for linear probing
            while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || strcmp((char*)s + h->keys[i] + 8, key) != 0)) {
                i = (i + inc) & mask;
                if (i == last) return h->n_buckets;
            }
            return __ac_iseither(h->flags, i)? h->n_buckets : i;
        } else return 0;
    }
    /* returns no characters processed. */
    unsigned encode_2bit(const uint8_t* seq)
    {
        size_t i = 0;
        uint8_t *p = dna, *rc = &dna[readlength-1];
        unsigned c, shft = 0;
        *rc = *p = '\0';
        cerr << seq << endl; //
        while ((c = *seq++) != '\0') {
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            *p |= c << shft; // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            *rc |= c << (shft - 6);
            shft = (++i & 0x3) << 1;
            if (shft == 0) *--rc = *++p = '\0';
        }
        return i;
    }
    uint8_t* encode_phredqual(const uint8_t* qual, const uint8_t* seq, uint8_t* b)
    {
        unsigned i = 0;
        uint8_t* p = b - 4;
        while ((*++b = *qual++) != '\0') {
            *b -= phred_offset;
            assert(*b < 40);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' 0x30 is set and phred bits originally there shifted upwards
            *b |= isb6(c) ? (c << 5) : ((*b & 0x30) << 2) | 0x30;
        }
        *p = (i & 0xff); i >>= 8; // insert phred sum
        *--p = (i & 0xff); i >>= 8;
        *--p = (i & 0xff); i >>= 8;
        *--p = (i & 0xff);
        return b;
    }
    void decode_phredqual(uint8_t* b, unsigned len)
    {
        uint8_t *d = dna, *q = dna + len + 1;

        for(unsigned i = 0; i != len; ++i, ++b, ++d, ++q) {
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
    }
    khash_t(UQCT) *H;
    size_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, keymax;
    static uint8_t *s;
    uint8_t* dna;
};

uint8_t *fq_hash::s;



#endif //RK_FQ_HASH_H

