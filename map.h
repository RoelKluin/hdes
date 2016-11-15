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

#ifndef RK_MAP_H
#define RK_MAP_H
#include "fa.h"

struct map_t {
    uint32_t p; // position and orientation
    uint32_t ho; // header offset
    uint32_t i; // location of k-mer within read
    uint32_t iend; // how many NTs were read in upseq
    uint32_t readlength;
    uint32_t mismatches; // or u8?
    uint32_t extension;
    uint32_t iteration;
    uint32_t nonk;
    uint8_t* s;
    char* upseq;
    //uint8_t* mism, unkey;
};


struct mapstat_t {
	unsigned keyct, read_too_long;
};

#endif // RK_MAP_H


