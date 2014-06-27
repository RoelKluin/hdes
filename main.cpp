/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
// make && zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 100000 | valgrind ./uqct

#include <stdint.h>
#include <assert.h>

#include <iostream>
#include <algorithm>

#include "khash.h"
#include "kstring.h"

#define MAX_NAME_ETC    1024

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

static unsigned phredtoe10[] = { 1, 1, 1, 1, 2, 3, 3, 5, 6, 7, 10,
    12, 15, 19, 25, 31, 39, 50, 63, 79, 100,
    125, 158, 199, 251, 316, 398, 501, 630, 794, 1000,
    1258, 1584, 1995, 2511, 3162, 3981, 5011, 6309, 7943, 10000 };

using namespace std;

#define decr4chtou32(b, i) ({   \
    i = *--b; i <<= 8;          \
    i |= *--b; i <<= 8;         \
    i |= *--b; i <<= 8;         \
    i | *--b;                   \
})

struct uniqct
{
    uniqct(size_t rl) : is_err(0), l(0), m(1ul << 23),
        fq_ent_max((rl << 1) + MAX_NAME_ETC), phred_offset(33), keymax(1)
    {
        H = kh_init(UQCT);
        s = (uint8_t*)malloc(m);
        *s = '\0';
        kroundup32(fq_ent_max);
    }

    void put(const uint8_t* seq, const uint8_t* qual, const uint8_t* name)
    {
        if_ever (is_err) return;
        if (l + fq_ent_max >= m) {
            m <<= 1;
            s = (uint8_t*)realloc(s, m);
        }
        uint8_t *p, *b = s + l;
        khiter_t k = str_kh_get(H, (const char*)seq);
        size_t i;
        if (k == kh_end(H)) { /* new key */
            for (p = b + 8; b != p; ++b) *b = '\0'; /* keycount/length */

            // TODO: form of twobit with N spec
            while ((*b = *seq) != '\0') ++b, ++seq;
            i = b - p;
            ++b;
            *--p = (i & 0xff); i >>= 8; /* length */
            *--p = (i & 0xff); i >>= 8;
            *--p = (i & 0xff);
//cerr << (unsigned)*p << "\t" << (unsigned)p[1] << "\t" << (unsigned)p[2] << "\t" << endl;
            i = p - s - 5; /* key offset, (denotes last value) */
            k = kh_put(UQCT, H, i, &is_err);
            if_ever ((is_err = (is_err < 0))) return;

        } else {
            p = s + kh_key(H, k);

            /* increment keycount */
            i = ++(*p);
            i |= (*++p += !i) << 8;
            i |= (*++p += !i) << 16;
            i |= (*++p += !i) << 24;
//cerr << "keycount:" << i << endl;
            if (i > keymax) keymax = i;
            i = kh_val(H, k); /* former value offset */
        }
        // 4 bytes reserved for phred sum
        b += 3;
        p = b;
        /* put former value offset here */
        for (unsigned j = 0; j != 4; ++j) {
            *++b = i & 0xff;
            i >>= 8;
        }
        // TODO: better qual/name storage
        i = 0;
        while ((*++b = *qual) != '\0') {
            assert(*qual - phred_offset < 40);
            i += phredtoe10[*qual - phred_offset];
            ++qual;
        }
        while ((*++b = *name) != '\0') ++name;
        l = ++b - s;
        kh_val(H, k) = p - s - 3;
        *p = (i & 0xff); i >>= 8; // insert phred sum
        *--p = (i & 0xff); i >>= 8;
        *--p = (i & 0xff); i >>= 8;
        *--p = (i & 0xff);

    }

    void dump() { /* must always be called to clean up */
        if (not is_err) {
            khiter_t *ks, *kp;
            size_t *vp, *vs;
            s[l] = '\0'; // end last value

            kp = ks = (khiter_t*)malloc(kh_size(H) * sizeof(khiter_t));
            for (khiter_t k = kh_begin(H); k != kh_end(H); ++k)
                if (kh_exist(H, k)) { *kp = k; ++kp; }
            ++keymax;
            kroundup32(keymax);
            cerr << "== Max number of identical sequences near " << keymax << endl;

            vs = (size_t*)malloc(keymax * sizeof(size_t));

            sort(ks, kp, *this); /* by keycount */

            //at = (uint64_t)((s[at+4] << 32)|(s[at+3] << 24)|(s[at+2] << 16)|(s[at+1] << 8)|s[at]);
            while (kp != ks) { /* loop over keys */
                size_t sq = kh_key(H, *--kp);
                size_t i = kh_val(H, *kp);
                uint8_t* p = s + sq + 5;
                unsigned len = *p; len <<= 8;
                len |= *++p; len <<= 8;
                len |= *++p;
                assert(len < (1 << 10)); //

/*first = sq, last = i; cerr << s + sq + 8 << "\ti:" << sq << "\t" <<i << endl; */

                vp = vs;
                do { // loop over values belonging to key
                    *vp = i;
                    ++vp;
                    uint8_t* b = s + i + 8;
                    i = decr4chtou32(b, i);

/* cerr << "i:" << i << endl; assert(i >= sq); */
                } while (i != sq);

                sort(vs, vp, *this); // by qualsum
                while (vp-- != vs) {
                    cout << (char*)s + *vp + 9 + len << "\n";
                    cout << (char*)s + sq + 8 << "\n+\n";
                    cout << (char*)s + *vp + 8 << "\n";
                }
            }
            free(vs);
            free(ks);
        }
	free(s);
        kh_destroy(UQCT, H);
    }
    inline bool operator() (khiter_t a, khiter_t b) /* needed by sort. No need to check length, no 2 keys the same.*/
    {
        uint8_t *sa = s + kh_key(H, a) + 8;
        uint8_t *sb = s + kh_key(H, b) + 8;
        while (*sa == *sb) ++sa, ++sb;
        return *sa > *sb;
    }

    inline bool operator() (size_t a, size_t b) /* needed by sort */
    {
/*assert(a >= first); assert(b >= first); assert(a <= last); assert(b <= last);*/
        uint8_t *sa = s + a, *sb = s + b;
        if (*sa != *sb) return *sa > *sb;
        if (*++sa != *++sb) return *sa > *sb;
        if (*++sa != *++sb) return *sa > *sb;
        return *++sa > *++sb;
    }
    int is_err;

private:
#define my_str_hash_equal(a, b) (strcmp((char*)s + a + 8, (char*)s + b + 8) == 0)
#define my_str_hash_func(key) __ac_X31_hash_string((char*)s + key + 8)
    KHASH_INIT2(UQCT, kh_inline, uint32_t, uint32_t, 1, my_str_hash_func, my_str_hash_equal)
    kh_inline khint_t str_kh_get(const kh_UQCT_t *h, const char* key)
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
//KHASH_MAP_INIT_STR(UQCT, uint32_t)
    khash_t(UQCT) *H;
    size_t l, m, fq_ent_max;
    unsigned phred_offset, keymax;
    uint8_t *s;
/*size_t first, last;*/
};


int main(int argc, const char* argv[])
{
    struct uniqct uqct(51);

    string nm, sq, tmp, ql;
    while ((not uqct.is_err) && getline(cin, nm) && getline(cin, sq) && getline(cin, tmp) && getline(cin, ql))
        uqct.put((const uint8_t*)sq.c_str(), (const uint8_t*)ql.c_str(), (const uint8_t*)nm.c_str());

    uqct.dump();

    if (uqct.is_err) {
        cerr << "An error occurred\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
