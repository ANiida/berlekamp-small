#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#include "global.h"
#include "struct.h"

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(array[0]))

//#define O  8  // 素数の和(3*5=15 の位数を持つ)
#define MAX  2  // 素数表の先頭から何個素数を足すか
#define NN   8  // 置換配列の次元

#define STR_LENGTH 128
#define PASSWORD_LENGTH 256

#define ROTL32(X, B) rotl32((X), (B))
static inline uint32_t
rotl32(const uint32_t x, const int b)
{
    return (x << b) | (x >> (32 - b));
}

#define ROTL64(X, B) rotl64((X), (B))
static inline uint64_t
rotl64(const uint64_t x, const int b)
{
    return (x << b) | (x >> (64 - b));
}

unsigned long xor128(void)
{
    static unsigned long x = 123456789, y = 362436069,
        z = 521288629, w = 88675123;
    unsigned int a = 0;
    unsigned long t;

    a = random();
    t = x ^ (a << 11);
    a = y;
    y = z;
    z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

int mlt(int x, int y)
{
    if (x == 0 || y == 0)
        return 0;

    return ((x + y - 2) % (G_M - 1)) + 1;
}

int mltn(int n, int x)
{
    int ret = 1;
    while (n > 0) {
        if (n & 1)
            ret = mlt(ret, x); // n の最下位bitが 1 ならば x^(2^i) をかける
        x = mlt(x, x);
        n >>= 1; // n を1bit 左にずらす
    }
    return ret;
}
