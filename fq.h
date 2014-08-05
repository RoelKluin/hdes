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

#ifndef RK_FQ_H
#define RK_FQ_H
#include "seq.h"
#include "gz.h"

typedef struct seqb2_t {
        uint8_t *s;
        uint32_t* lookup;
	uint64_t mode, l, m;
        uint32_t nr, key_ct, readlimit;
        uint16_t readlength;
        uint8_t phred_offset; //XXX:
        struct gzfh_t fh[4]; /* reader and writer */
} seqb2;

int init_fq(seqb2_t *fq);
int fq_b2(seqb2_t *fq);
void free_fq(seqb2_t *fq);
int fq_print(seqb2_t *fq);
#endif //FQ_ARR_H
