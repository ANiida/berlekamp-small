#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "8192.h"

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 8 // 符号長
#define M 8 //有限体の元の数
#define K (4) //符号の次元
#define E (3)    //拡大体のビット数
#define DEG K*E //(K * E)
#define T (K / 2) //エラーの数

