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

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

// get bit for alpha char, regardless of case [command line options]
#define amopt(r)        (1ul << ((r+7) & 0x1f))

/* the following should be enough for 32 bit int */

static inline uint8_t *
sprints(uint8_t *out, uint8_t *s)
{
    while (*s) { *++out = *s; ++s; }
    return out;
}


static inline unsigned
sprintu(uint8_t *out, register uint64_t u, uint8_t del)
{
    register uint8_t *s = out + 21;
    *s = del;

    while (u) {
        register unsigned t = u % 10;
        u /= 10;
        *--s = t + '0';
    }
    for (;*s != del; ++u, ++s) *++out = *s;
    return ++u;
}

static inline unsigned
sprint0x(uint8_t *out, register uint64_t u, uint8_t del, uint8_t len)
{
    register unsigned i;
    *++out = '0'; *++out = 'x';
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
