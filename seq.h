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

#ifndef RK_SEQ_H
#define RK_SEQ_H
#include "gz.h"
#include "b6.h"
#include "util.h"

#define SEQ_OFFSET_BYTES 2
#define NS_OFFSET_BYTES 3
#define BUF_OFFSET_BYTES 5
#define SEQ_OFFSET_MAX ((1ul << (SEQ_OFFSET_BYTES << 3)) - 1ul)

#define INIT_BUFSIZEBIT 23
#define INIT_BUFSIZE (1ul << INIT_BUFSIZEBIT)

typedef struct seqb2_t {
        uint8_t * s;
        uint32_t* lookup;
	uint64_t mode, s_l;
        uint32_t nr, key_ct, maxreads;
        uint16_t readlength, blocksize;
        uint8_t phred_offset, s_m; //XXX: bitfields?
        struct gzfh_t fh[5]; /* reader and writer */
} seqb2;

unsigned get_fastx_type(char*, const unsigned);

#endif //RK_SEQ_H
