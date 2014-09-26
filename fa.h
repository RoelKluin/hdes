/*
 *    Description:
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel Kluin,
 */

#ifndef RK_FA_H
#define RK_FA_H
#include "seq.h"
#include "gz.h"
#define BUF_STACK (1 << 22)

typedef struct index_action {
    const char* search;
    const char* replace;
    int (*index) (seqb2_t*, void*, int (*) (void*));
} index_action_t;

int fa_index(seqb2_t *seq);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
