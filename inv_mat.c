//
// (over finite field) Gauss-Jordan法による逆行列
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"
#include "struct.h"
#include "gf.h"

int mlt(int x, int y);    /// in chash.c

int Inv(unsigned short b)
{
    if (b == 0)
        return 0;

    for (int i = 0; i < G_N; i++) {
        if (gf[mlt(i, b)] == 1)
            return i;
    }
    return -1;
}

MTX mulmat(MTX A, MTX B, int flg)  //// flg ????
{
    MTX tmp = {0};
    /*
      switch (flg) { ????? 1, 2, 3 で何が違う ????
      ループ回数も、回し方も違う？？？
      1:  tmp.x[j][i] ^= A.x[i][k] & B.x[j][k];
      2:  tmp.x[j][i] ^= A.x[i][k] & B.x[j][k];
      3:  tmp.x[i][j] ^= gf[mlt(fg[A.x[i][k]], fg[B.x[k][j]])];
    */
    switch (flg) {
    case 1:
        for (int i = 0; i < G_DEG; i++)    {
            for (int j = 0; j < G_N; j++)      {
                for (int k = 0; k < G_DEG; k++)        {
                    tmp.x[j][i] ^= A.x[i][k] & B.x[j][k];
                }
            }
        }
        break;
    case 2:
        for (int i = 0; i < G_AY; i++)    {
            for (int j = 0; j < G_N; j++)      {
                for (int k = 0; k < G_AY; k++)        {
                    tmp.x[j][i] ^= A.x[i][k] & B.x[j][k];
                }
            }
        }
        break;
    case 3:
        for (int i = 0; i < G_DEG; i++)    {
            for (int j = 0; j < G_DEG; j++)      {
                for (int k = 0; k < G_DEG; k++)        {
                    tmp.x[i][j] ^= gf[mlt(fg[A.x[i][k]], fg[B.x[k][j]])];
                }
            }
        }
        break;
    }

    return tmp;
}
