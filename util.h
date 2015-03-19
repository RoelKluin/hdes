/*
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel KLuin, 
 */

#ifndef RK_UTIL_H
#define RK_UTIL_H
#include <stdint.h>
#include <assert.h>


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

// the comma before the ## is deleted in absense of an argument
#define WARN(msg, ...) EPR("Warning: " msg " at %s:%u", ##__VA_ARGS__, __FILE__, __LINE__)

#define ASSERT(cond, action, ...) \
if (!(cond)) { \
    WARN("assertion '" #cond "' failed " __VA_ARGS__);\
    action;\
}

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

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
#define _buf_init(buf, sz) _buf_init_err(buf, sz, return -ENOMEM)

#define _buf_init_arr(buf, sz) ({\
    buf##_m = sz;\
    /*fprintf(stderr, #buf " malloc, %lu\n", sizeof(*(buf)) << buf##_m);fflush(NULL);*/\
    decltype(buf) __t = (decltype(buf))malloc(sizeof(*(buf)) << buf##_m);\
    if_ever (__t == NULL) return -ENOMEM;\
    __t;\
})

#define _buf_grow_err_print(buf, step, error_action) \
if (buf##_l + step >= (1ul << buf##_m)) {\
    fprintf(stderr, #buf " realloc, %lu\n", sizeof(*(buf)) << (buf##_m + 1));fflush(NULL);\
    decltype(buf) __t = (decltype(buf))realloc(buf, sizeof(*(buf)) << ++(buf##_m));\
    if_ever (__t == NULL) error_action;\
    buf = __t;\
}

#define _buf_grow_err(buf, step, error_action) \
if (buf##_l + step >= (1ul << buf##_m)) {\
    /*fprintf(stderr, #buf " realloc, %lu\n", sizeof(*(buf)) << (buf##_m + 1));fflush(NULL);*/\
    decltype(buf) __t = (decltype(buf))realloc(buf, sizeof(*(buf)) << ++(buf##_m));\
    if_ever (__t == NULL) error_action;\
    buf = __t;\
}
#define _buf_grow(buf, step) _buf_grow_err(buf, step, return -ENOMEM)
#define _buf_grow_add_err(buf, step, add, error_action) ({\
    _buf_grow_err(buf, step, error_action);\
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
sprintu(uint8_t *out, register uint64_t u, uint8_t del)
{
    register unsigned len = 0;
    register uint8_t *s = out + 21;

    do {
        ++len;
        register unsigned t = u % 10;
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
sprint0x(uint8_t *out, register uint64_t u, uint8_t del, uint8_t len)
{
    register unsigned i;
    *out = '0'; *++out = 'x';
    register uint8_t *s = out + len + 1;
    *s = del;

    // hex: zero padded.
    for (i = 0; i != len; ++i) {
        register int t = u & 0xf;
        *--s = t + (t >= 10 ? 'a' - 10 : '0');
        u >>= 4;
    }
    return len + 3;
}




#define DEBUG 1
#endif // RK_UTIL_H
