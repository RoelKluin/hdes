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

// maximum length of the read name
#define SEQ_MAX_NAME_ETC    (1u << 14)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ_SHFT ((KEY_LENGTH * 2) - 1)
#define KEYNT_BUFSZ (1ul << KEYNT_BUFSZ_SHFT)
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1ul)

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_MASK ((1u << (KEY_WIDTH << 1)) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)
#define KEYNT_TRUNC_UPPER (~HALF_KEYNT_MASK & KEYNT_TRUNC_MASK)


#define SNDX_TRUNC_MASK (KEYNT_TRUNC_MASK | KEYNT_BUFSZ)
#define HALF_SNDX_MASK (HALF_KEYNT_MASK | KEYNT_STRAND)
#define SNDX_TRUNC_UPPER (~HALF_SNDX_MASK & SNDX_TRUNC_MASK)

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
        struct gzfh_t fh[4]; /* reader and writer */
} seqb2;

unsigned get_fastx_type(char* f, const unsigned last_fq, const unsigned fhsz);


#endif //RK_SEQ_H
