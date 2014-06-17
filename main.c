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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "khash.h"
KHASH_MAP_INIT_INT64(READ, uint64_t)

#define kh_count(n, k, h, key, ret, on_err)     \
    do {                                        \
        k = kh_get(n, h, key);                  \
        if (k == kh_end(h)) {                   \
            k = kh_put(n, h, key, &ret);        \
            if (ret == -1) on_err;              \
            kh_val(h, k) = 1;                   \
        } else {                                \
            kh_val(h, k)++;                     \
        }                                       \
    } while (0);

int main(int argc, const char* argv[])
{
    int i, ret;
    khiter_t k;
    uint64_t key;

    khash_t(READ) *h = kh_init(READ);

    for (i = 1; i != argc; ++i) {
        key = atol(argv[i]);
        fprintf(stderr, "orig:%s\t%lu\n", argv[i], key);
        fflush(NULL);
        kh_count(READ, k, h, key, ret, return EXIT_FAILURE);
    }

    for (k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k))
            fprintf(stderr, "%lu\t%lu\n", kh_key(h, k), kh_val(h, k));
    }

    kh_destroy(READ, h);
    return EXIT_SUCCESS;
}
