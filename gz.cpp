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
            int i;
            fprintf(stderr, "ERROR:%s:%d:", gzerror(fh->io, &i), c);
            fprintf(stderr, "%d\n", i);
            return c;
        }
        //EPQ(dbg > 5, "==%d bytes written", c);
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
            int i;
            fprintf(stderr, "ERROR:%s:%d:", gzerror(fh->io, &i), c);
            fprintf(stderr, "%d\n", i);

            return c;
        }
        //fprintf(stderr, "==%d bytes read\n", c);
        l -= c, s += c;
    }

    return gzeof(fh->io) ? -EIO : 0;
}

int
b2_write(const gzfh_t *fh, const char *s, uint64_t l)
{
    while (l) {
        int c = write(fileno(fh->fp), s, l > INT_MAX ? INT_MAX : l);
        if (c < 0) {
            EPR("error while writing");
            return c;
        }
        //EPR("==%d bytes written", c);
        l -= c, s += c;
    }
    return ferror(fh->fp) ? -3 : 0;
}

static int
b2_read(const gzfh_t *fh, char *s, uint64_t l)
{
    while (l) {
        int c = read(fileno(fh->fp), s, l > INT_MAX ? INT_MAX : l);
        if (c < 0) {
            EPR("error while reading");
            return c;
        }
        //EPR("==%d bytes read", c);
        l -= c, s += c;
    }
    return ferror(fh->fp) ? -3 : 0;
}

void set_readfunc(struct gzfh_t* fhin, void** g, int (**gc)(void*),
        int (**ungc)(int, void*))
{
    if (fhin->io == NULL) {
        *g = fhin->fp;
        *gc = (int (*)(void*))&fgetc;
        *ungc = (int (*)(int, void*))&ungetc;
    } else {
        *g = fhin->io;
        *gc = (int (*)(void*))&gzgetc;
        *ungc = (int (*)(int, void*))&gzungetc;
    }
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
rgzopen(struct gzfh_t* fh)
{
    /* wb1: 1 is for fast compression */
    fh->io = gzdopen(fileno(fh->fp), fh->write ? "wb1" : "rb");

    if (fh->io == Z_NULL) {
        fputs("main: gzdopen failed\n", stderr);
        return -EACCES;
    }

    if (gzbuffer(fh->io, fh->blocksize << 20) == -1) {
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

static int
fn_convert(struct gzfh_t* fhio, const char* search, const char* replace)
{
    char* f = strstr(fhio->name, search);
    ASSERT(f != NULL, return -EFAULT,
            "Could not find extension %s (and convert to %s)", search, replace);
    strncpy(f, replace, strlen(replace));
    return 0;
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
        if (fh->fp != stdout && fh->fp != stdin) {
            EPR("closing fh->fp");
            fclose(fh->fp);
        }
    } else {
        EPQ(dbg > 6, "already seems closed: %s", fh->name ? fh->name : "<unknown>");
    }
    fh->fp = NULL;
    return ret;
}

// a force of 2 forces read, return 0 for reading, for writing 1, else error.
int
set_io_fh(struct gzfh_t* fh, int force)
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
        fh->read = NULL;
        if (t && t[3] == '\0') {/* '.gz' at _end_ */
            fh->write = b2gz_write;
            ret = rgzopen(fh);
        } else {
	    fh->write = b2_write;
	}
        return ret;
    } else {
        fprintf(stderr, "== preparing to read %s...\n", fh->name);
        if (fh->fp == NULL) {
            fprintf(stderr, "cannot read from %s\n", fh->name);
            return -ENOENT;
        }
	fh->write = NULL;
        if (t && t[3] == '\0') {
            fh->read = b2gz_read;
            return rgzopen(fh);
        } else {
            fh->read = b2_read;
	}
    }
    return 0;
}

// open file if it exists. Zero means ok, test fhio->fp to see whether file existed.
int
reopen(struct gzfh_t* fhio, const char* search, const char* replace)
{
    int res = rclose(fhio);
    if (res < 0)
        return res;
    _ACTION(fn_convert(fhio, search, replace), "")

    fhio->fp = fopen(fhio->name, "r");
    res = 0;
err:
    return res;
}



