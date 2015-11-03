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
    gzFile io;
    int (*read) (const gzfh_t *fh, char*, size_t);
    int (*write) (const gzfh_t *fh, const char*, size_t);
    int fd;
    uint16_t blocksize;
} gzfh;
int b2_write(const gzfh_t*, const char*, size_t);
int set_stdio_fh(struct gzfh_t* fh, uint64_t* mode);
int set_io_fh(struct gzfh_t* fh, int force);
int rclose(gzfh_t *fh);
int reopen(struct gzfh_t*, const char*, const char*);
void set_readfunc(struct gzfh_t*, void**, int (**)(void*));
#endif //RK_GZ_H
