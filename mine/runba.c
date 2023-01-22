#include <stdint.h>
//#include "struct.h"
#include "8192.h"
#include "ranbada.c"

unsigned short gf_mul(unsigned short in0, unsigned short in1)
{
    int i;

    uint32_t tmp;
    uint32_t t0;
    uint32_t t1;
    uint32_t t;

    t0 = in0;
    t1 = in1;

    tmp = t0 * (t1 & 1);

    for (i = 1; i < 12; i++)
        tmp ^= (t0 * (t1 & (1 << i)));

    t = tmp & 0x7FC000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    t = tmp & 0x3000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    return tmp & ((1 << 12) - 1);
}

/* input: in0, in1 in GF((2^m)^t)*/
/* output: out = in0*in1 */
void GF_mul(unsigned short *out, unsigned short *in0, unsigned short *in1)
{
    int i, j;

    unsigned short prod[K * 2 - 1] = {0};

    for (i = 0; i < K * 2 - 1; i++)
        prod[i] = 0;

    for (i = 0; i < K; i++)
    {
        for (j = 0; j < K; j++)
            prod[i + j] ^= gf[mlt(fg[in0[i]], fg[in1[j]])];
    }
    //

    for (i = (K - 1) * 2; i >= K; i--)
    {
        if (K == 512)
        {
            // GF(2^512) from sage
            prod[i - K + 8] ^= prod[i];
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 256)
        {
            // GF(2^256) from sage
            prod[i - K + 10] ^= prod[i];
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 128)
        {
            // 128
            prod[i - K + 7] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 1] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 32)
        {
            // 32
            prod[i - K + 15] ^= prod[i];
            prod[i - K + 9] ^= prod[i];
            prod[i - K + 7] ^= prod[i];
            prod[i - K + 4] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 16)
        {
            // 16
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 8)
        {
            // 8
            prod[i - K + 4] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 4)
        {
            // 4
            prod[i - K + 1] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
    }

    for (i = 0; i < K; i++)
        out[i] = prod[i];
}
