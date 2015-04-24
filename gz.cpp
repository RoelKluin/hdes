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
#include <limits.h>
#include "gz.h"
static inline int
b2gz_write(const gzfh_t* fh, const char* s, size_t l)
{
    while (l) {
        // last argument is unsigned
        int c = gzwrite(fh->io, s, l > UINT_MAX ? UINT_MAX : l);
        if (c < 0) {
            fprintf(stderr, "%s\n", gzerror(fh->io, &c));
            return c;
        }
        //fprintf(stderr, "==%d bytes written\n", c);
        l -= c, s += c;
    }

    return gzeof(fh->io) ? -EIO : 0;
}

static inline int
b2gz_read(const gzfh_t* fh, char* s, size_t l)
{
    while (l) {
        // last argument is unsigned
        int c = gzread(fh->io, s, l > UINT_MAX ? UINT_MAX : l);
        if (c < 0) {
            fprintf(stderr, "%s\n", gzerror(fh->io, &c));
            return c;
        }
        //fprintf(stderr, "==%d bytes read\n", c);
        l -= c, s += c;
    }

    return gzeof(fh->io) ? -EIO : 0;
}

static int
rgzopen(struct gzfh_t* fh, uint32_t blocksize)
{
    int fd = fileno(fh->fp);
    fh->fd = dup(fd);
    /* wb1: 1 is for fast compression */
    fh->io = gzdopen(fd, fh->write ? "wb1" : "rb");

    if (fh->io == Z_NULL) {
        fputs("main: gzdopen failed\n", stderr);
        return -EACCES;
    }

    if (gzbuffer(fh->io, blocksize << 20) == -1) {
        fputs("main: gzdopen failed to set blocksize\n", stderr);
        return -EINVAL;
    }
    return 0;
}

int
set_stdio_fh(struct gzfh_t* fh, uint64_t* mode)
{
    const char* t = strstr(fh->name, fh->write ? "stdout" : "stdin");
    if (!t) return 0;

    if (*mode & amopt(t[3])) {
        fprintf(stderr, "%s was already assigned\n", t);
        return -EINVAL;
    }
    fh->fp = fh->write ? stdout : stdin;
    *mode |= amopt(t[3]);
    return 1;
}

/* a force of 2 forces read, return 0 for reading 1, writing, else error. */
int
set_io_fh(struct gzfh_t* fh, uint32_t blocksize, int force)
{
    if (fh->fp == NULL)
        fh->fp = fopen(fh->name, "r"); // test whether file exists
    const char* t = strstr(fh->name,".gz");

    if ((fh->write && force != 2) || fh->fp == NULL) { /* write without gzip by default */
        fprintf(stderr, "== preparing to write %s...\n", fh->name);
        if (fh->fp) { /* file already exists, overwrite? */
            fclose(fh->fp);
            if (!force) {
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
            return rgzopen(fh, blocksize) + 1;
        }
    } else {
        fprintf(stderr, "== preparing to read %s...\n", fh->name);
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot read from %s\n", fh->name);
            return -ENOENT;
        }
        if (t && t[3] == '\0') {
            fh->read = b2gz_read;
            fh->close = gzclose_r;
            return rgzopen(fh, blocksize);
        }
    }
    return 0;
}



