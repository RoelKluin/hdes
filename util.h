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

#ifndef RK_UTIL_H
#define RK_UTIL_H
#include <stdint.h>
#include <assert.h>
#include <execinfo.h>
#include <stdlib.h>
#include <signal.h>
#include "b6.h"

#define DEBUG 1

#define C const

#ifndef kroundup
// g++ may complain if you place shift after test.
#define kroundup(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4,\
                            (x)|=(x)>>(sizeof(x) > 1 ? 8 : 0),\
                            (x)|=(x)>>(sizeof(x) > 2 ? 16 : 0),\
                            (x)|=(x)>>(sizeof(x) > 4 ? 32 : 0),\
                            (x)|=(x)>>(sizeof(x) > 8 ? 64 : 0),\
                            ++(x))
#endif

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#define STR(s) #s

#define EPQ(d, msg, ...) if (d) { fprintf (stderr, msg "\n", ##__VA_ARGS__); }
#define EPQ0(d, ...) if (d) { fprintf (stderr, __VA_ARGS__); }
#define EPR(msg, ...) do { fprintf (stderr, msg "\n", ##__VA_ARGS__); fflush(stderr); } while(0)
#define EPR0(...) do { fprintf (stderr, __VA_ARGS__); fflush(stderr); } while(0)
#define OPR(msg, ...) fprintf (stdout, msg "\n", ##__VA_ARGS__)
#define OPR0(...) fprintf (stdout, __VA_ARGS__)

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(x)    if (unlikely(x))
#define expect(x)       if (likely(x))
#define while_ever(x)    while (unlikely(x))

// the comma before the ## is deleted in absense of an argument
#define WARN(msg, ...) EPR("Warning: " msg " at %s:%u", ##__VA_ARGS__, __FILE__, __LINE__)
#define QARN(d, msg, ...) EPQ(d, "Warning: " msg " at %s:%u", ##__VA_ARGS__, __FILE__, __LINE__)

/* Note: The "MS" section flags are to remove duplicates.  */
#define DEFINE_GDB_PY_SCRIPT(script_name) \
 asm("\
.pushsection \".debug_gdb_scripts\", \"MS\",@progbits,1\n\
.byte 1 /* Python */\n\
.asciz \"" gdb_py/frame_filter.py "\"\n\
.popsection \n\
");


static void handler(int sig) {
  void *array[50];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 50);

  // print out all the frames to stderr
  fprintf(stderr, "\nError: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

// continue debugging after assertion failure.
// http://stackoverflow.com/questions/1721543/continue-to-debug-after-failed-assertion-on-linux#1721575
#if DEBUG > 0
#define NB(cond, ...) \
do {\
    if_ever (!(cond)) { \
        WARN("assertion '" #cond "' failed " __VA_ARGS__);\
        raise(SIGTRAP);\
    }\
} while(0)
#else
#define NB(cond, ...) (cond)
#endif

#define ASSERT(cond, action, ...) \
do {\
    if_ever (!(cond)) { \
        WARN("assertion '" #cond "' failed " __VA_ARGS__);\
        action;\
    }\
} while(0)

#define _ACTION0(fun, msg, ...)\
do {\
    res = fun;\
    NB(res >= 0, msg "(%u) %s:%u.", ##__VA_ARGS__, res, __FILE__, __LINE__);\
    EPQ(msg[0] != '\0', msg ".", ##__VA_ARGS__);\
} while(0)

#define _EVAL(fun)\
    NB((res = (fun)) >= 0, "(%u) at %s:%u.", res, __FILE__, __LINE__);


#define _ACTION(fun, msg, ...)\
do {\
    res = fun;\
    NB(res >= 0, msg "(%u) at %s:%u.", ##__VA_ARGS__, res, __FILE__, __LINE__);\
    EPQ(msg[0] != '\0', msg "..\tdone", ##__VA_ARGS__);\
} while(0)

#define packed_struct struct __attribute__ ((__packed__))

// get bit for alpha char, regardless of case [command line options]
#define amopt(r)        (1ul << ((r+7) & 0x1f))
#define buf_init(buf, sz) ({\
    buf##_m = sz;\
    buf##_l = 0;\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << buf##_m);\
    if_ever (__t == NULL)\
        raise(SIGTERM);\
    __t;\
})

#define buf_init_arr(buf, sz) ({\
    buf##_m = sz;\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << buf##_m);\
    if_ever (__t == NULL)\
        raise(SIGTERM);\
    __t;\
})

#define buf_grow_shift(buf, step, shft) \
do {\
    decltype(buf##_l) __t = ((buf##_l + step) >> shft);\
    if (__t >= (1ul << buf##_m)) {\
        kroundup(__t);\
        buf##_m = __builtin_ctz(__t) + 1;\
        buf = (decltype(buf))realloc(buf, sizeof(*(buf)) << buf##_m);\
        if_ever (buf == NULL)\
            raise(SIGTERM);\
    }\
} while(0)

#define buf_grow(buf, step) buf_grow_shift(buf, step, 0)

#define buf_grow_add(buf, step, push) ({\
    buf_grow_shift(buf, step, 0);\
    buf[buf##_l++] = push;\
})

#define buf_grow_struct_add(buf, ...) ({\
    buf_grow_shift(buf, 1ul, 0);\
    buf[buf##_l++] = { __VA_ARGS__ };\
})



#define buf_free(buf) \
do {\
    if (buf != NULL) {\
        free(buf);\
        /*buf = NULL;*/\
    }\
} while (0)
/* copy buffer and return pointer to one past it
 * s must be zero terminated, and out have sufficient size always.
 */
static inline uint8_t *
sprints(uint8_t *out, uint8_t *s)
{
    for (; *s; ++out, ++s) *out = *s;
    return out;
}

/*
 * almost equivalent to `sprints(out, "%lu%c", u, del)'
 * the following should be enough for 64 bit int
 */
static inline unsigned
sprintu(uint8_t *out, uint64_t u, uint8_t del)
{
    unsigned len = 0;
    uint8_t *s = out + 21;

    do {
        ++len;
        unsigned t = u % 10;
        u /= 10;
        *--s = t + '0';
    } while (u);
    while (u++ != len) *out++ = *s++;
    *out = del;
    return len + 1;
}

/*
 * almost equivalent to `sprints(out, "0x%016lx%c", u, del)' if len == 16.
 * again should be enough for 64 bit int
 */
static inline unsigned
sprint0x(uint8_t *out, uint64_t u, uint8_t del, uint8_t len)
{
    unsigned i;
    *out = '0'; *++out = 'x';
    uint8_t *s = out + len + 1;
    *s = del;

    // hex: zero padded.
    for (i = 0; i != len; ++i) {
        int t = u & 0xf;
        *--s = t + (t >= 10 ? 'a' - 10 : '0');
        u >>= 4;
    }
    return len + 3;
}

static unsigned next_pow2(unsigned x)
{
    x -= 1;
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    return x + 1;
}

static void
print_dna(uint32_t dna, const char sep = '\n', unsigned len = KEY_WIDTH)
{
    for (unsigned t = len; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc(sep, stderr);
}
static void
print_2dna(uint32_t dna, uint32_t dna2, unsigned len = KEY_WIDTH)
{
    print_dna(dna, '|', len);
    print_dna(dna2, '\n', len);
}
static void
print_seq(keyseq_t* seq, unsigned len = KEY_WIDTH)
{
    print_dna(seq->dna, '|', len);
    print_dna(seq->rc, '\n', len);
}

static void
print_posseq(uint8_t const*const s, uint32_t p, unsigned len = KEY_WIDTH)
{
    keyseq_t seq {.p = p = _b2pos_of(p)};
    build_seq(s, seq, p);
    EPR0("%u\t", p>>1);
    print_dna(seq.dna, '\n', len);
}

static void
print_ndx(uint32_t dna)
{
    uint32_t rc = dna & KEYNT_TRUNC_UPPER;
    dna ^= rc ^ (rc << 1) ^ KEYNT_STRAND;
    rc = revcmp(dna);
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('|', stderr);
    for (unsigned t = KEY_WIDTH; t--; rc >>= 2)
        fputc(b6((rc & 3) << 1), stderr);
    fputc('\n', stderr);
}


// to get array size of array member of t
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

// Internal Macros
#define HEX__(n) 0x##n##LU
#define B8__(x) ((x&0x0000000FLU)?1:0) \
  +((x&0x000000F0LU)?2:0) \
  +((x&0x00000F00LU)?4:0) \
  +((x&0x0000F000LU)?8:0) \
  +((x&0x000F0000LU)?16:0) \
  +((x&0x00F00000LU)?32:0) \
  +((x&0x0F000000LU)?64:0) \
  +((x&0xF0000000LU)?128:0)

// User-visible Macros
#define B8(d) ((unsigned char)B8__(HEX__(d)))
#define B16(dmsb,dlsb) (((unsigned short)B8(dmsb)<<8) + B8(dlsb))
#define B32(dmsb,db2,db3,dlsb) \
  (((unsigned long)B8(dmsb)<<24) \
  + ((unsigned long)B8(db2)<<16) \
  + ((unsigned long)B8(db3)<<8) \
  + B8(dlsb))

#endif // RK_UTIL_H
