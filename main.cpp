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
#include <iostream>
#include <algorithm>

#include "khash.h"
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expected(T)     if (likely(T))

using namespace std;

struct uniqct {
    uniqct(char* buf, size_t bl) : b(buf), bstart(b), bend(b + bl), buflen(bl)
    {
        h = kh_init(UQCT);
        is_err = (h == NULL);
    }
    bool ok() { return not is_err; }
    void put(const char* key)
    {
        if_ever (is_err) return;
        khiter_t k = kh_get(UQCT, h, key);
        if (k == kh_end(h)) {
            size_t l = strlen(key) + 1; /* include '\0' */
            is_err = (b + l >= bend);
            expected (ok()) {
                strncpy(b, key, l);
                k = kh_put(UQCT, h, b, &is_err); //
                is_err = (is_err < 0);
                b += l;
            }
            if_ever (is_err) return;
            kh_val(h, k) = 1;
        }
        kh_val(h, k)++;
    }
    void dump() {
        unsigned l;
        khiter_t *ks;
        if_ever (is_err) goto err;
        ks = (khiter_t*)malloc(kh_size(h) * sizeof(khiter_t));
        is_err = (ks == NULL);
        if_ever (is_err) goto err;
        l = 0;
        for (khiter_t k = kh_begin(h); k != kh_end(h); ++k)
            if (kh_exist(h, k)) ks[l++] = k;

        std::sort(ks, ks + l, *this);

        while (l--)
            cout << kh_val(h, ks[l]) << "\t" << kh_key(h, ks[l]) << endl;

        free(ks);
err:	kh_destroy(UQCT, h);

    }
    bool operator() (khiter_t a, khiter_t b) /* needed by sort */
    {
        return kh_val(h, a) != kh_val(h, b) ?
            kh_val(h, a) > kh_val(h, b) : strcmp(kh_key(h, a), kh_key(h, b)) > 0;
    }
private:
    KHASH_INIT2(UQCT, kh_inline, const char*, uint32_t, 1, kh_str_hash_func, kh_str_hash_equal)
    char *b, *bstart, *bend;
    khash_t(UQCT) *h;
    uint64_t buflen;
    int is_err;
};

int main(int argc, const char* argv[])
{
    unsigned buflen = 1ul << 31;
    char* buf = (char*)malloc(buflen * sizeof(char));
    struct uniqct uqct(buf, buflen);

    string line;
    while (uqct.ok() && cin && getline(cin, line))
        uqct.put(line.c_str());

    uqct.dump();
    free(buf);

    return uqct.ok() ? EXIT_SUCCESS : EXIT_FAILURE;
}
