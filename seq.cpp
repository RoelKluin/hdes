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
 * Return type of extension. 0:
 */
unsigned
get_fastx_type(char* f, const unsigned ref_i, const unsigned fhsz)
{
    unsigned i = 0, c = strlen(f) - 1;
    f += c; // parse extension from right

    if (*f == 'z' && (((c -= 3) < 3) || *--f != 'g' || *--f != '.')) //.gz
        return fhsz;

    if (*--f == 't') { // .txt(.gz)?
        if (((c -= 4) < 0) || *--f != 'x' || *--f != 't') return fhsz + 1;
    } else {
        if (*f == 'a') i = ref_i;
        else if (*f != 'q') return fhsz + 2;

        if (*--f == 't') { // .fast[aq](.gz)??
            if (((c -= 3) < 3) || *--f != 's' || *--f != 'a') return fhsz + 3;
            --f;
        }
        if (*f != 'f') return fhsz + 4;
    }
    return *--f == '.' ? i : fhsz + 5;
}


