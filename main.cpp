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
#include <stdint.h>
#include <iostream>
#include <algorithm>

#include "khash.h"

using namespace std;
KHASH_MAP_INIT_INT(RD, uint32_t)

struct uniqct {
    uniqct() : ks(NULL), maxlen(0), at(0) { h = kh_init(RD); }
    int put(uint32_t key)
    {
        khiter_t k = kh_get(RD, h, key);
        if (k == kh_end(h)) {
            int ret = 0;
            k = kh_put(RD, h, key, &ret);
            if ((ret < 0) || ((at >= maxlen) && ks_realloc(maxlen << 1))) {
                at = 0;
                dump();
                return -1;
            }
            ks[at++] = k;
            kh_val(h, k) = 0;
        } else {
            kh_val(h, k)++;
        }
        return 0;
    }
    int ks_realloc(int ml) {
        if (ml == 0) ml = 1024;
        khiter_t* t = (khiter_t*)realloc(ks, ml * sizeof(khiter_t));
        if (t == NULL) return -1;
        ks = t;
        maxlen = ml;
        return 0;
    }

    void sort() { if (ks) std::sort(ks, ks + at, *this); }
    void dump() {
        while (at--)
            cout << kh_val(h, ks[at]) << "\t" << kh_key(h, ks[at]) << endl;
        free(ks);
        kh_destroy(RD, h);
    }
    bool operator() (khiter_t a, khiter_t b)
    {
        return kh_val(h, a) != kh_val(h, b) ?
            kh_val(h, a) < kh_val(h, b) : kh_key(h, a) < kh_key(h, b);
    }
private:
    khash_t(RD) *h;
    khiter_t* ks;
    int maxlen, at;
};

int main(int argc, const char* argv[])
{
    struct uniqct uqt;
    
    string line;
    while (cin) {
        getline(cin, line);
        if (uqt.put(atoi(line.c_str())) < 0)
            return EXIT_FAILURE;
    }
    //for (int i = 1; i != argc; ++i)

    //std::sort(ks, ks + kc.size(), h);
    uqt.sort();
    uqt.dump();

    return EXIT_SUCCESS;
}
