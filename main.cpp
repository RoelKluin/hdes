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
#include "kstring.h"

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expected(T)     if (likely(T))

using namespace std;

struct uniqct {
    uniqct(char* b, kstring_t* s, size_t l, size_t n) :
        Bstart(b), B(b), Bend(b + l),
        Sstart(s), S(s), Send(s + n)
    {
        H = kh_init(UQCT);
        is_err = (H == NULL);
    }
    void put(const char* key, const char* val)
    {
        if_ever (is_err) return;
        khiter_t k = kh_get(UQCT, H, key);
        if (k == kh_end(H)) {
            size_t l = strlen(key) + 1; /* include '\0' */
            is_err = (B + l >= Bend);
            expected (ok()) {
                strncpy(B, key, l);
                k = kh_put(UQCT, H, B, &is_err); //
                is_err = (is_err < 0);
                B += l;
            }
            if_ever (is_err) return;
            S->l = S->m = 0;
            S->s = NULL;
            kputs(val, S);
            kh_val(H, k) = S++ - Sstart;
        } else {
            kputs(val, &Sstart[kh_val(H, k)]);
        }
    }
    void dump() {
        unsigned l;
        khiter_t *ks;
        if_ever (is_err) goto err;
        ks = (khiter_t*)malloc(kh_size(H) * sizeof(khiter_t));
        is_err = (ks == NULL);
        if_ever (is_err) goto err;
        l = 0;
        for (khiter_t k = kh_begin(H); k != kh_end(H); ++k)
            if (kh_exist(H, k)) ks[l++] = k;

        std::sort(ks, ks + l, *this);

        while (l--)
            cout << kh_key(H, ks[l]) << "\n" << Sstart[kh_val(H, ks[l])].s;

        free(ks);
err:	while (S-- != Sstart) free(S->s);
        kh_destroy(UQCT, H);

    }
    bool operator() (khiter_t a, khiter_t b) /* needed by sort */
    {
        size_t ma = Sstart[kh_val(H, a)].m;
        size_t mb = Sstart[kh_val(H, b)].m;
        return ma != mb ? ma > mb : strcmp(kh_key(H, a), kh_key(H, b)) > 0;
    }
    bool ok() { return not is_err; }
private:
    KHASH_INIT2(UQCT, kh_inline, const char*, uint32_t, 1, kh_str_hash_func, kh_str_hash_equal)
    char *Bstart, *B, *Bend;
    kstring_t *Sstart, *S, *Send;
    khash_t(UQCT) *H;
    int is_err;
};

int main(int argc, const char* argv[])
{
    unsigned buflen = 1ul << 30, splen = 1ul << 28; // 268,435,456 reads, 
    char* buf = (char*)malloc(buflen * sizeof(char));
    kstring_t* sp = (kstring_t*)malloc(buflen * sizeof(kstring_t));
    struct uniqct uqct(buf, sp, buflen, splen);

    string nm, sq, np, ql;
    while (uqct.ok() && getline(cin, nm) && getline(cin, sq) && getline(cin, np) && getline(cin, ql)) {
        nm += "\t" + ql + "\n";
        uqct.put(sq.c_str(), nm.c_str());
    }

    uqct.dump();
    free(buf);
    free(sp);

    return uqct.ok() ? EXIT_SUCCESS : EXIT_FAILURE;
}
