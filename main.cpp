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
    i = (uint8_t)*--b; i <<= 8;          \
    i |= (uint8_t)*--b; i <<= 8;         \
    i |= (uint8_t)*--b; i <<= 8;         \
    i | (uint8_t)*--b;                   \
})

struct uniqct
{
    uniqct(size_t rl) : is_err(0), l(0),
        m(1ul << 23), fq_ent_max((rl << 1) + MAX_NAME_ETC), keymax(1)
    {
        H = kh_init(UQCT);
        s = (char*)malloc(m);
        *s = '\0';
        kroundup32(fq_ent_max);
    }

    void put(const char* key, const char* val)
    {
        if_ever (is_err) return;
        if (l + fq_ent_max >= m) {
            m <<= 1;
            s = (char*)realloc(s, m);
        }
        char *p, *b = s + l;
cerr << key << endl << flush;
        khiter_t k = str_kh_get(H, key);
        size_t i;
        if (k == kh_end(H)) { /* new key */
            p = b;
            for (i = 0; i != 5; ++i, ++b) *b = '\0'; /* keycount */
            for (i = 0; i != 3; ++i, ++b) *b = '\0'; /* seq length */

            // TODO: form of twobit with N spec
            while ((*b++ = *key) != '\0') ++key;
            i = p - s; /* key offset, (denotes last value) */
            k = kh_put(UQCT, H, i, &is_err);
            is_err = (is_err < 0);
            if_ever (is_err) return;

cerr << "key at: " << i << endl << "kh_size: " << kh_size(H) << endl << flush;

        } else {
            p = s + kh_key(H, k);

cerr << "already observed key at " << kh_key(H, k) << ":\n" << p + 8 << "\n";

            /* increment keycount */
            i = (size_t)++(*p);
            i |= (size_t)(*++p += !i) << 8;
            i |= (size_t)(*++p += !i) << 16;
            keymax |= ((size_t)(*++p += !i) << 24) | i;
            i = kh_val(H, k); /* former value offset */
        }

        p = b;
        // 4 bytes reserved for phred sum
        *b = '\0'; *++b = '\0'; *++b = '\0'; *++b = '\0'; ++b;
        /* put former value offset here */
        for (unsigned j = 0; j != 4; ++j, ++b) {
            *b = i & 0xff;
cerr << i << ":" << (i & 0xff) << endl;
            i >>= 8;
        }
        // TODO: better qual/name storage
        while ((*b = *val) != '\0') ++b, ++val;
        kh_val(H, k) = p - s;
cerr << p + 8 << endl << "val at: " << p - s << endl << flush;
        l = ++b - s;
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

            sort(ks, kp, *this);

            //at = (uint64_t)((s[at+4] << 32)|(s[at+3] << 24)|(s[at+2] << 16)|(s[at+1] << 8)|s[at]);
            while (kp != ks) { /* loop over keys */
                size_t end = kh_key(H, *--kp);
                size_t i = kh_val(H, *kp);
                cout << s + end + 8 << ":key\n";
                cout << s + i + 8 << ":val\n";

                vp = vs;
                do { /* loop over values belonging to key */
cerr << "offsetted key/val at: " << i + 8 << endl;
                    *vp = i;
                    ++vp;
                    char* b = s + i + 8;
                    i = decr4chtou32(b, i);

cerr << "got key/val: " << i << "/" << end << endl;

if (i > (1 << 12)) { cerr << "ERROR\n"; break; }
if (vp - vs != 0 && vp[-1] == i) { cerr << "ERROR2\n"; break; }
if ((vp - vs > 1) && vp[-2] == i) { cerr << "ERROR3\n"; break; }
if (vp - vs > keymax) { cerr << "ERROR4\n"; break; }

                } while (i != end);

if (i > (1 << 12)) break;
//if (vp - vs != 1 && vp[-2] == i) break;
//if ((vp - vs > 2) && vp[-3] == i) { cerr << "ERROR5\n"; break; }
//if (vp - vs >= keymax) { cerr << "ERROR4\n"; break; }

cerr << "sorting " << vp - vs << " keys ...\n";
                sort(vs, vp, *this);

                while (vp != vs) cout << s + *--vp + 8 << "\n";
            }
            free(vs);
            free(ks);
        }
	free(s);
        kh_destroy(UQCT, H);
    }
    inline bool operator() (khiter_t a, khiter_t b) /* needed by sort */
    {
        char *sa = s + kh_key(H, a) + 8;
        char *sb = s + kh_key(H, b) + 8;
        while (*sa == *sb) { ++sa; ++sb; }
        return *sa > *sb;
    }

    inline bool operator() (size_t a, size_t b) /* needed by sort */
    {
        return strcmp(s + a, s + b);
    }
    int is_err;

private:
#define my_str_hash_equal(a, b) (strcmp(s + a + 8, s + b + 8) == 0)
#define my_str_hash_func(key) __ac_X31_hash_string(s + key + 8)
    KHASH_INIT2(UQCT, kh_inline, uint32_t, uint32_t, 1, my_str_hash_func, my_str_hash_equal)
    kh_inline khint_t str_kh_get(const kh_UQCT_t *h, const char* key)
    {
        if (h->n_buckets) {
            khint_t inc, k, i, last, mask;
            mask = h->n_buckets - 1;
            k = __ac_X31_hash_string(key); i = k & mask;
            inc = __ac_inc(k, mask); last = i; // inc==1 for linear probing
            while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || strcmp(s + h->keys[i] + 8, key) != 0)) {
                i = (i + inc) & mask;
                if (i == last) return h->n_buckets;
            }
            return __ac_iseither(h->flags, i)? h->n_buckets : i;
        } else return 0;
    }
//KHASH_MAP_INIT_STR(UQCT, uint32_t)
    khash_t(UQCT) *H;
    size_t l, m, fq_ent_max, keymax;
    char *s;
};


int main(int argc, const char* argv[])
{
    bool order_by_qual = true;
    // bool discard_comment = true;
    // bool entrylength_varies_widely = false;
    struct uniqct uqct(51);

    string nm, sq, tmp, ql;
    while ((not uqct.is_err) && getline(cin, nm) && getline(cin, sq) && getline(cin, tmp) && getline(cin, ql)) {
        tmp = order_by_qual ? ql + "\t" + nm : nm + "\t" + ql;
        uqct.put(sq.c_str(), tmp.c_str());
    }

    uqct.dump();

    if (uqct.is_err) {
        cerr << "An error occurred\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
