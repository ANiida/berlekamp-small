//#include "struct.h"
#include "8192.h"
#include "ranbada.c"

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
    int i, n; //, pid;
    vec s = {0}, u = {0}, r = {0};
    vec v = {0}; //, ff=o2v(f);
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    // int m = E;
    //  m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = (v);
    // for (i = 0; i < K / 2; i++)
    r = s;
    n = deg((f));

    if (vLT(f).n == 0 && vLT(f).a == 1)
    {
        printf("f==0\n");
        exit(1);
    }
    if (n == 0)
        return -1;

    i = 0;

    // r(x)^{q^i} square pow mod
    for (i = 0; i < K / 2; i++)
    {
        printf(":i=%d", i);
        // irreducible over GH(8192) 2^13
        r = vpp(r, f);
        // if(r.x[0]==65535)
        // return -1;
        u = vadd(r, (s));
        u = vgcd(f, u);

        if (deg(u) > 0 || vLT(u).a == 0)
        {
            // flg[i]= -1;
            printf("ae\n");
            return -1;
        }
    }

    return 0;
}

// #define NN 16
vec renritu(MTX a)
{
    unsigned short p, d;
    int i, j, k;
    vec v = {0};

    for (i = 0; i < K; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K + 1); j++)
        {
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (j = 0; j < K; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K + 1); k++)
                {
                    a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
                }
            }
        }
    }
    for (i = 0; i < K; i++)
    {
        if (a.x[i][i] != 1)
            // exit(1);
            for (j = 0; j < K + 1; j++)
                printf("%d,", a.x[i][j]);
        printf("\n");
    }
    printf("\n");

    for (i = 0; i < K; i++)
    {
        v.x[i] = a.x[i][K];
        // v.x[128]=1;
        printf(" x%d = %d\n", i, v.x[i]);
    }

    return v;
}

