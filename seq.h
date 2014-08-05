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


// maximum length of the read name
#define FQ_MAX_NAME_ETC    (1u << 14)

// KEY_LENGTH must be odd - or 2nd bit of central Nt is not always flipped in its
// complement - the alternative is the twisted halfdev conversion, but this is cheaper

#define KEY_LENGTH 15                       // <= this many Nts are used as key
#define KEY_CUTOFF 4                        // cut off, to get random keys
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF) // entire key maximized on (after conversion)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ (1u << (KEY_LENGTH * 2 - 1))
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1u)

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_TOP ((KEY_WIDTH - 1) * 2)
#define KEYNT_MASK ((1ul << KEY_WIDTH * 2) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)

#define SEQ_OFFSET_BYTES 2
#define BUF_OFFSET_BYTES 4
#define SEQ_OFFSET_MAX ((1u << (SEQ_OFFSET_BYTES << 3)) - 1u)

#endif //RK_SEQ_H

