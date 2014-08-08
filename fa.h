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

int init_fa(seqb2_t *fa);
int fa_ndx(seqb2_t *fa);
void free_fa(seqb2_t *fa);
int fa_print(seqb2_t *fa);
#endif // RK_FA_H
