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
#include "b6.h"

#define DEBUG 1

static const unsigned long dbgndx = ~0ul;// 0x1d03c;//0x1c27012;
static unsigned dbgkct = -3u;
static const uint32_t dbgndxkct = -3u; //1099511627775; //0x2028;
static const char* dbgrn = "HWI-ST745_0097:7:1101:7550:1094#0/1";
static const char* dbgchr = "B";//GL000226.1";//GL000207.1T";//GL000197.1";//GL000239.1";//GL000231.1";
static const unsigned long dbgtsoffs = 5000;//11;//212469;//2687783;//4995737;//3964212;//1835528;//3888783;//1955549;//2811688;//~0ul;
static const uint32_t dbgpos = 0; //1099511627775; //0x2028;
static int dbg = 3;

#define C const

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#define STR(s) #s

#define EPQ(d, msg, ...) if (d) { fprintf (stderr, msg "\n", ##__VA_ARGS__); fflush(NULL); }
#define EPQ0(d, ...) if (d) { fprintf (stderr, __VA_ARGS__); fflush(NULL); }
#define EPR(msg, ...) fprintf (stderr, msg "\n", ##__VA_ARGS__)
#define EPR0(...) fprintf (stderr, __VA_ARGS__)
#define OPR(msg, ...) fprintf (stdout, msg "\n", ##__VA_ARGS__)
#define OPR0(...) fprintf (stdout, __VA_ARGS__)

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

// the comma before the ## is deleted in absense of an argument
#define WARN(msg, ...) EPR("Warning: " msg " at %s:%u", ##__VA_ARGS__, __FILE__, __LINE__)
#define QARN(d, msg, ...) EPQ(d, "Warning: " msg " at %s:%u", ##__VA_ARGS__, __FILE__, __LINE__)

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

#define ASSERT(cond, action, ...) \
if_ever (!(cond)) { \
    WARN("assertion '" #cond "' failed " __VA_ARGS__);\
    action;\
}

#define _ACTION0(fun, msg, ...)\
    res = fun;\
    ASSERT(res >= 0, goto err, msg "(%u) %s:%u.", ##__VA_ARGS__, res, __FILE__, __LINE__);\
    EPQ(msg[0] != '\0', msg ".", ##__VA_ARGS__);

#define _EVAL(fun)\
    ASSERT((res = (fun)) >= 0, goto err, "(%u) at %s:%u.", res, __FILE__, __LINE__);


#define _ACTION(fun, msg, ...)\
    res = fun;\
    ASSERT(res >= 0, goto err, msg "(%u) at %s:%u.", ##__VA_ARGS__, res, __FILE__, __LINE__);\
    EPQ(msg[0] != '\0', msg "..\tdone", ##__VA_ARGS__);

#define packed_struct struct __attribute__ ((__packed__))

// get bit for alpha char, regardless of case [command line options]
#define amopt(r)        (1ul << ((r+7) & 0x1f))
#define _buf_init_err(buf, sz, error_action) ({\
    buf##_m = sz;\
    buf##_l = 0;\
    /*fprintf(stderr, #buf " malloc, %lu\n", sizeof(*(buf)) << buf##_m);fflush(NULL);*/\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << buf##_m);\
    if_ever (__t == NULL) error_action;\
    __t;\
})
#define _buf_init_err_m(buf, m, error_action) ({\
    buf##_l = 0;\
    /*fprintf(stderr, #buf " malloc, %lu\n", sizeof(*(buf)) << buf##_m);fflush(NULL);*/\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << m);\
    if_ever (__t == NULL) error_action;\
    __t;\
})
#define _buf_init(buf, sz) _buf_init_err(buf, sz, return -ENOMEM)

#define _buf_init_arr_err(buf, sz, error_action) ({\
    buf##_m = sz;\
    /*fprintf(stderr, #buf " malloc, %lu\n", sizeof(*(buf)) << buf##_m);fflush(NULL);*/\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << buf##_m);\
    if_ever (__t == NULL) error_action;\
    __t;\
})

#define _buf_grow_err_m(buf, step, m, shft, error_action) \
if (((buf##_l + step) >> shft) >= (1ul << m)) {\
    EPQ(dbg < 0, #buf " realloc, %lu", sizeof(*(buf)) << (buf##_m + 1));\
    decltype(buf) __t = (decltype(buf))realloc(buf, sizeof(*(buf)) << ++m);\
    if_ever (__t == NULL) error_action;\
    buf = __t;\
}
#define _buf_grow(buf, step, shft) _buf_grow_err_m(buf, step, buf##_m, shft, return -ENOMEM)
#define _buf_growm(buf, step, m, shft) _buf_grow_err_m(buf, step, m, shft, return -ENOMEM)
#define _buf_grow0(buf, step) _buf_grow_err_m(buf, step, buf##_m, 0, return -ENOMEM)
#define _buf_grow0m(buf, step, m) _buf_grow_err_m(buf, step, m, 0, return -ENOMEM)

#define _buf_grow_err(buf, step, shft, err) _buf_grow_err_m(buf, step, buf##_m, shft, err)
#define _buf_grow_add_err(buf, step, shft, add, error_action) ({\
    _buf_grow_err(buf, step, shft, error_action);\
    buf[buf##_l++] = add;\
})



#define _buf_grow2(buf, step, s) \
if (buf##_l + step >= (1ul << buf##_m)) {\
    /*fprintf(stderr, #buf " realloc, %lu\n", sizeof((*buf)) << (buf##_m + 1));fflush(NULL);*/\
    s = (decltype(buf))realloc(buf, sizeof((*buf)) << ++(buf##_m));\
    if_ever (s == NULL) return -ENOMEM;\
    buf = s;\
    s += buf##_l;\
}

#define _buf_free(buf) \
if (buf != NULL) {\
    free(buf);\
    /*buf = NULL;*/\
}
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

static int
print_dna(uint64_t dna, bool d = true, const char sep = '\n', unsigned len = KEY_WIDTH)
{
    if (d) {
        for (unsigned t = len; t--; dna >>= 2)
            fputc(b6((dna & 3) << 1), stderr);
        fputc(sep, stderr);
    }
    return -1;
}
static int
print_2dna(uint64_t dna, uint64_t dna2, bool d = true, unsigned len = KEY_WIDTH)
{
    if (d) {
        print_dna(dna, true, '|', len);
        print_dna(dna2, true, '\n', len);
    }
    return -1;
}

static int
print_ndx(uint64_t dna, bool d = true)
{
    if (d == false) return -1;
    uint64_t rc = dna & KEYNT_TRUNC_UPPER;
    dna ^= rc ^ (rc << 1) ^ KEYNT_STRAND;
    rc = revcmp(dna);
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('|', stderr);
    for (unsigned t = KEY_WIDTH; t--; rc >>= 2)
        fputc(b6((rc & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
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

#define DEBUG 1
#endif // RK_UTIL_H
