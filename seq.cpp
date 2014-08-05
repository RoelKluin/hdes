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
#include <stdlib.h> // malloc()
#include <string.h> // memcpy()
#include <errno.h> // ENOMEM
#include "seq.h"
/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
int
init_seq(seqb2_t *seq)
{
    seq->m = INIT_BUFSIZE;
    seq->s = (uint8_t*)malloc(INIT_BUFSIZE);
    if (seq->s == NULL) return -ENOMEM;
    *seq->s = '\0';

    size_t i, v = 1, l = sizeof(*seq->lookup) * KEYNT_BUFSZ;
    seq->lookup = (uint32_t*)malloc(l);
    if (seq->lookup == NULL) {
        free(seq->s);
        return -ENOMEM;
    }
    for (i = 0; i != l; i += sizeof(*seq->lookup))
        memcpy(((char*)seq->lookup) + i, &v, sizeof(*seq->lookup));
    return 0;
}

void
free_seq(seqb2_t *seq)
{
    if (seq->lookup != NULL) { free(seq->lookup); seq->lookup = NULL; }
    if (seq->s != NULL) { free(seq->s); seq->s = NULL; }
}



