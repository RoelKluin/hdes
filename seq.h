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

// KEY_WIDTH must be odd - or 2nd bit of central Nt is not always flipped in its
// complement - the alternative is the twisted halfdev conversion, but this is cheaper

#define KEY_LENGTH 16                       // <= this many Nts are used as key
#define KEY_CUTOFF 5                        // cut off, to get random keys, odd!
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF) // entire key maximized on (after conversion)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ (1ul << (KEY_LENGTH * 2))
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1u)

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_TOP ((KEY_WIDTH - 1) * 2)
#define KEYNT_MASK ((1ul << KEY_WIDTH * 2) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)

#define SEQ_OFFSET_BYTES 2
#define BUF_OFFSET_BYTES 4
#define SEQ_OFFSET_MAX ((1u << (SEQ_OFFSET_BYTES << 3)) - 1u)

#define INIT_BUFSIZE (1u << 23)

typedef struct seqb2_t {
        uint8_t *s;
        uint32_t* lookup;
	uint64_t mode, l, m;
        uint32_t nr, key_ct, readlimit;
        uint16_t readlength, blocksize;
        uint8_t phred_offset; //XXX:
        struct gzfh_t fh[4]; /* reader and writer */
} seqb2;

unsigned get_fastx_type(char* f, const unsigned last_fq, const unsigned fhsz);


#endif //RK_SEQ_H
