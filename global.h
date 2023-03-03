#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 符号のパラーメータの指定。通常[N,K,T]として、
// Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
// を表す。ここではDは符号長にしている。

#if 1
/// これはちゃんと動く。
# define G_N    8192       // 符号長 (== 256)
# define G_M    8192       // 有限体の元の数 (M <= N)
# define G_K    128        // 符号の次元
# define G_E    13         // 拡大体のビット数
#else
/// ulimit -s unlimited しても、いきなりセグフォ。
# define G_N   8192       // 符号長
# define G_M   8192       // 有限体の元の数
# define G_K    256       // 符号の次元
# define G_E     13       // 拡大体のビット数
#endif

#define G_DEG     (G_K * G_E)   /// 次数 (degree)
#define G_K2      (G_K / 4)     /// 削除した関数で使われていた
#define G_T       (G_K / 2)     // エラーの数
///    for (j = 1; j < G_T * 2 + 1; j++) {  BK.c:1256

#define G_AY     ((G_K / 2 + 1) * G_E)      /// 元は lu.c にあった "AY"
#define G_K21     (G_K / 2 + 1)
///#define G_K21E ((G_K / 2 + 1) * G_E)
