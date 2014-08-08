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
#include <errno.h>
#include <string.h> //strstr()
#include "gz.h"

static inline int
b2gz_write(gzfh_t *fh, char *start, int len)
{
    if ((len = gzwrite(fh->io, start, len)) < 0) {
        fprintf(stderr, "%s\n", gzerror(fh->io, &len));
        return len;
    }
    if (DEBUG)
        fprintf(stderr, "%d bytes written\n", len);

    return gzeof(fh->io) ? -EIO : len;
}

static int
rgzopen(struct gzfh_t* fh, uint64_t blocksize)
{
    int fd = fileno(fh->fp);
    fh->fd = dup(fd);
    /* wb1: 1 is for fast compression */
    fh->io = gzdopen(fd, fh->write ? "wb1" : "rb");

    if (fh->io == Z_NULL) {
        fputs("main: gzdopen failed\n", stderr);
        return -EACCES;
    }

    if (gzbuffer(fh->io, blocksize) == -1) {
        fputs("main: gzdopen failed to set blocksize\n", stderr);
        return -EINVAL;
    }
    return 0;
}

int
set_io_fh(struct gzfh_t* fh, uint64_t* mode, uint64_t blocksize)
{
    const char* t = strstr(fh->name, fh->write ? "stdout" : "stdin");
    if (t) {
        if (*mode & amopt(t[3])) {
            fprintf(stderr, "%s was already assigned\n", t);
            return -EINVAL;
        }
        fh->fp = fh->write ? stdout : stdin;
        *mode |= amopt(t[3]);
        return 0;
    }
    fh->fp = fopen(fh->name, "r"); // test whether file exists
    t = strstr(fh->name,".gz"); // TODO: read magic gzip number in file

    if (fh->write) { /* write without gzip by default */
        if (fh->fp) { /* file already exists, overwrite? */
            fclose(fh->fp);
            if ((*mode & amopt('f')) == 0) {
                fprintf(stderr, "use -f to force write %s\n", fh->name);
                return -EEXIST;
            }
        }
        fh->fp = fopen(fh->name, "w");
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot write to %s\n", fh->name);
            return -EPERM;
        }
        if (t && t[3] == '\0') {/* '.gz' at _end_ */
            fh->write = b2gz_write;
            fh->close = gzclose_w;
            return rgzopen(fh, blocksize);
        }
    } else {
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot read from %s\n", fh->name);
            return -ENOENT;
        }
        if (t && t[3] == '\0') {
            fh->close = gzclose_r;
            return rgzopen(fh, blocksize);
        }
    }
    return 0;
}



