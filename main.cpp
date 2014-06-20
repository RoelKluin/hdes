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
            is_err = (B + l >= Bend) | (S + 1 >= Send);
            if_ever (is_err) return;

            strncpy(B, key, l);
            k = kh_put(UQCT, H, B, &is_err); //
            is_err = (is_err < 0);
            if_ever (is_err) return;

            B += l;
            S->l = S->m = 0;
            S->s = NULL;
            kputs(val, S);
            kh_val(H, k) = S++ - Sstart;

        } else {
            kputs(val, &Sstart[kh_val(H, k)]);
        }
    }

    void dump() {
        khiter_t *ks, *kp;
        if (is_err) goto err;

        ks = (khiter_t*)malloc(kh_size(H) * sizeof(khiter_t));
        is_err = (ks == NULL);
        if (is_err) goto err;

        kp = ks;
        for (khiter_t k = kh_begin(H); k != kh_end(H); ++k)
            if (kh_exist(H, k)) *kp++ = k;

        sort(ks, kp, *this);

        while (kp-- != ks)
            cout << kh_key(H, *kp) << "\n" << Sstart[kh_val(H, *kp)].s;

        free(ks);
err:	while (S-- != Sstart) free(S->s);
        if (H) kh_destroy(UQCT, H);
        H = NULL;
    }

    bool operator() (khiter_t a, khiter_t b) /* needed by sort */
    {
        size_t xa, xb;
        xa = Sstart[kh_val(H, a)].m; /* sort by number of entries */
        xb = Sstart[kh_val(H, b)].m;
        if (xa != xb) return xa > xb; /* maxes compare slightly faster since they are powers of 2 */

        xa = Sstart[kh_val(H, a)].l;
        xb = Sstart[kh_val(H, b)].l;
        return (xa == xb) ? strcmp(kh_key(H, a), kh_key(H, b)) > 0 : xa > xb;
    }
    int is_err;

private:
    KHASH_MAP_INIT_STR(UQCT, uint32_t)
    char *Bstart, *B, *Bend;
    kstring_t *Sstart, *S, *Send;
    khash_t(UQCT) *H;
};

int main(int argc, const char* argv[])
{
    size_t buflen = 1ul << 28, splen = 1ul << 26; // 67,108,864 reads
    char* buf = (char*)malloc(buflen * sizeof(char));
    kstring_t* sp = (kstring_t*)malloc(splen * sizeof(kstring_t));

    bool order_by_qual = true;
    // bool discard_comment = true;
    // bool entrylength_varies_widely = false;
    struct uniqct uqct(buf, sp, buflen, splen);

    string nm, sq, tmp, ql;
    while ((not uqct.is_err) && getline(cin, nm) && getline(cin, sq) && getline(cin, tmp) && getline(cin, ql)) {
        tmp = (order_by_qual ? ql + "\t" + nm : nm + "\t" + ql) + "\n";
        uqct.put(sq.c_str(), tmp.c_str());
    }

    uqct.dump();
    free(buf);
    free(sp);

    if (uqct.is_err) {
        cerr << "An error occurred\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
