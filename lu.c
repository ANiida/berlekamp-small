// 学校で習った知識を使ったら問題解決できました。大学行っといてよかったｗ

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"
#include "struct.h"

MTX inv_S = {0};
MTX norm_S = {0};    /// MTX S = {0};

int is_reg(MTX cc, MTX *R)
{
    unsigned int flg_count = 0, count;  /// flg_count は bool ではない！
    unsigned char    cl[G_DEG][G_DEG];
    unsigned char     b[G_DEG][G_DEG] = {0};
    unsigned char inv_a[G_DEG][G_DEG] = {0}; //ここに逆行列が入る

    while (flg_count != G_DEG) {     //// ????
        time_t t;

        flg_count = 0;             //// 0 から？
        srand(clock() + time(&t));

        for (int i = 0; i < G_DEG; i++) {
            for (int j = 0; j < G_DEG; j++) {
                cl[i][j] = cc.x[i][j];
            }
        }

        // 単位行列を作る
        for (int i = 0; i < G_DEG; i++) {
            for (int j = 0; j < G_DEG; j++) {
                inv_a[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        //掃き出し法
        for (int i = 0; i < G_DEG; i++) {
            if (cc.x[i][i] == 0) {
                int j = i;

                while (cc.x[j][i] == 0 && j < G_DEG) {
                    j++;
                }

                if (j > G_DEG) {
                    printf("S is not reg in is_reg %d\n", j);
                        return -1;
                }

                for (int k = 0; k < G_DEG; k++) {
                    cc.x[i][k] ^= cc.x[j][k] % 2;
                    inv_a[i][k] ^= inv_a[j][k];
                }

                cc.x[i][i] = 1;
            }

            if (cc.x[i][i] == 1)  {
                for (int l = i + 1; l < G_DEG; l++) {
                    if (cc.x[l][i] == 1) {
                        for (int k = 0; k < G_DEG; k++) {
                            cc.x[l][k] ^= cc.x[i][k] % 2;
                            inv_a[l][k] ^= inv_a[i][k];
                        }
                    }
                }
            }
        }

        for (int i = 1; i < G_DEG; i++) {
            for (int k = 0; k < i; k++) {
                if (cc.x[k][i] == 1) {
                    for (int j = 0; j < G_DEG; j++) {
                        cc.x[k][j] ^= cc.x[i][j] % 2;
                        inv_a[k][j] ^= inv_a[i][j];
                    }
                }
            }
        }

        //検算
        for (int i = 0; i < G_DEG; i++) {
            for (int j = 0; j < G_DEG; j++)  {
                for (int k = 0; k < G_DEG; k++) {    ///// ここまで直した
                    b[i][j] ^= (cl[i][k] & inv_a[k][j]);
                }
            }
        }

        for (int i = 0; i < G_DEG; i++) {
            if (b[i][i] == 1) {
                flg_count++;
            }
        }

        count = 0;
        for (int i = 0; i < G_DEG; i++) {
            for (int j = 0; j < G_DEG; j++) {
                if (b[i][j] == 0 && i != j)
                    count++;
            }
        }
        printf("%d,%d %d,%d", flg_count, G_DEG, count, G_DEG * G_DEG - G_DEG);

        printf("S[K][K]=\n{\n");
        if (flg_count == G_DEG && count == (G_DEG * G_DEG - G_DEG)) {
            printf("inv_S[K][K]=\n{\n");
            for (int i = 0; i < G_DEG; i++) {
                printf("{");
                for (int j = 0; j < G_DEG; j++) {
                    R->x[i][j] = inv_a[i][j];
                    printf("%d,", inv_S.x[i][j]);
                }
                printf("},\n");
            }
            printf("};\n");
            return 0;
        }
        return -1;
    }
    return -1;
}

int mkRS(MTX cc, MTX *R)
{
    unsigned char b[G_AY][G_AY] = {0};
    unsigned int flg = 0;                 //// bool ではない
    unsigned int count;     
    unsigned char cl[G_AY][G_AY];
    unsigned char inv_a[G_AY][G_AY] = {0};    //ここに逆行列が入る

    while (flg != G_AY) {
        time_t t;

        flg = 0;
        count = 0;
        srand(clock() + time(&t));

        printf("end of g2\n");

        for (int i = 0; i < G_AY; i++) {
            for (int j = 0; j < G_AY; j++)  {
                cl[i][j] = cc.x[i][j];
            }
        }

        //単位行列を作る
        for (int i = 0; i < G_AY; i++)  {
            for (int j = 0; j < G_AY; j++)    {
                inv_a[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        //掃き出し法
        for (int i = 0; i < G_AY; i++)    {
            if (cc.x[i][i] == 0)      {
                int j = i;

                while (cc.x[j][i] == 0 && j < G_AY) {
                    j++;
                }

                if (j >= G_AY)  {
                    printf("baka in mkS %d\n", j);
                    return -1;
                }
                for (int k = 0; k < G_AY; k++) {
                    cc.x[i][k] ^= cc.x[j][k] % 2;
                    inv_a[i][k] ^= inv_a[j][k];
                }

                cc.x[i][i] = 1;
            }

            if (cc.x[i][i] == 1)            {
                for (int l = i + 1; l < G_AY; l++)  {
                    if (cc.x[l][i] == 1)      {
                        for (int k = 0; k < G_AY; k++) {
                            cc.x[l][k] ^= cc.x[i][k] % 2;
                            inv_a[l][k] ^= inv_a[i][k];
                        }
                    }
                }
            }
        }

        for (int i = 1; i < G_AY; i++) {
            for (int k = 0; k < i; k++) {
                if (cc.x[k][i] == 1) {
                    for (int j = 0; j < G_AY; j++) {
                        cc.x[k][j] ^=  cc.x[i][j] % 2;
                        inv_a[k][j] ^= inv_a[i][j];
                    }
                }
            }
        }

        //検算
        for (int i = 0; i < G_AY; i++) {
            for (int j = 0; j < G_AY; j++) {
                for (int k = 0; k < G_AY; k++) {
                    b[i][j] ^= (cl[i][k] & inv_a[k][j]);
                }
            }
        }

        for (int i = 0; i < G_AY; i++) {
            if (b[i][i] == 1) {
                flg++;
            }
        }
        count = 0;

        for (int i = 0; i < G_AY; i++) {
            for (int j = 0; j < G_AY; j++) {
                if (b[i][j] == 0 && i != j)
                    count++;
            }
        }

        printf("norm_S[K][K]=\n{\n");    //// printf("S[K][K]=\n{\n");
        if (flg == G_AY && count == (G_AY * G_AY - G_AY)) {
            for (int i = 0; i < G_AY; i++) {
                for (int j = 0; j < G_AY; j++) {
                    norm_S.x[i][j] = cl[i][j];
                    printf("%d,", norm_S.x[i][j]);
                }
                printf("},\n");
            }
            printf("};\n");

            printf("inv_S[K][K]=\n{\n");
            for (int i = 0; i < G_AY; i++) {
                for (int j = 0; j < G_AY; j++) {
                    R->x[i][j] = inv_a[i][j];
                    printf("%d,", inv_S.x[i][j]);    //// アヤシイ。
                }
            }
            printf("},\n");
            return 0;
        }
        return -1;
    }
    return -1;
}
