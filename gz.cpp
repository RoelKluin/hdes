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
        int c = gzwrite(fh->io, s, l > INT_MAX ? INT_MAX : l);
        if (c <= 0) {
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
        int c = gzread(fh->io, s, l > INT_MAX ? INT_MAX : l);
        if (c < 0) {
            fprintf(stderr, "%s\n", gzerror(fh->io, &c));
            return c;
        }
        //fprintf(stderr, "==%d bytes read\n", c);
        l -= c, s += c;
    }

    return gzeof(fh->io) ? -EIO : 0;
}

/* from: http://www.zlib.net/manual.html at gzdopen:
 * gzclose on the gzopen-returned gzFile will also close the file descriptor fd.
 *
 * If you want to keep fd open, use fd = dup(fd_keep); gz = gzdopen(fd, mode);.
 * The duplicated descriptor should be saved to avoid a leak, since gzdopen does
 * not close fd if it fails. If you are using fileno() to get the file descriptor
 * from a FILE *, then you will have to use dup() to avoid double-close()ing the
 * file descriptor. Both gzclose() and fclose() will close the associated file
 * descriptor, so they need to have different file descriptors.
 */
static int
rgzopen(struct gzfh_t* fh, uint32_t blocksize)
{
    /* wb1: 1 is for fast compression */
    fh->io = gzdopen(fileno(fh->fp), fh->write ? "wb1" : "rb");

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

int
rclose(gzfh_t* fh)
{
    int ret = 0;
    if (fh->io) {
        if (fh->read && fh->write) {
            ret = gzclose(fh->io);
        } else if (fh->read) {
            ret = gzclose_r(fh->io);
            fclose(fh->fp);
        } else if (fh->write) {
            ret = gzclose_w(fh->io);
            fclose(fh->fp);
        } else if (fh->fp) {
            fclose(fh->fp);
        }
        fh->read = NULL;
        fh->write = NULL;
        fh->io = NULL;
    } else if (fh->fp) {
        EPR("closing fh->fp");
        fclose(fh->fp);
    } else {
        EPR("already seems closed: %s", fh->name ? fh->name : "<unknown>");
    }
    fh->fp = NULL;
    return ret;
}

// a force of 2 forces read, return 0 for reading, for writing 1, else error.
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
                fh->fp = NULL;
                fprintf(stderr, "use -f to force write %s\n", fh->name);
                return -EEXIST;
            }
        }
        fh->fp = fopen(fh->name, "w");
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot write to %s\n", fh->name);
            return -EPERM;
        }
        int ret = 1;
        if (t && t[3] == '\0') {/* '.gz' at _end_ */
            fh->write = b2gz_write;
            ret = rgzopen(fh, blocksize);
        }
        return ret;
    } else {
        fprintf(stderr, "== preparing to read %s...\n", fh->name);
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot read from %s\n", fh->name);
            return -ENOENT;
        }
        if (t && t[3] == '\0') {
            fh->read = b2gz_read;
            return rgzopen(fh, blocksize);
        }
    }
    return 0;
}



