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

#ifndef RK_GZ_H
#define RK_GZ_H
#include <stdio.h>
#include <zlib.h>
#include "util.h"

typedef struct gzfh_t {
    char* name;
    FILE *fp;
    int fd;
    gzFile io;
    int (*write) (const gzfh_t *fh, const char*, size_t);
    int (*close) (gzFile);
    uint64_t blocksize;
} gzfh;

int set_stdio_fh(struct gzfh_t* fh, uint64_t* mode);
int set_io_fh(struct gzfh_t* fh, uint32_t blocksize, int force);

#endif //RK_GZ_H
