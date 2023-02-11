#include <stdlib.h>

#include "global.h"
//#include "8192.h"
//#include "4096.h"
//#include "512.h"
//#include "256.h"
//#include "16.h"

// static const unsigned short gf[8]={0,1,2,4,3,6,7,5};
// static const unsigned short fg[8]={0,1,2,4,3,7,5,6};

// nomal bases
static unsigned short gf[G_M] = { 0, 1, 2, 4, 8, 9, 11, 15, 7, 14, 5, 10, 13, 3, 6, 12 };
static unsigned short fg[G_M] = { 0, 1, 2, 13, 3, 10, 14, 8, 4, 5, 11, 6, 15, 12, 9, 7 };

// sage比較用
// static const unsigned short gf[16]={0,1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};
// static const unsigned short fg[16]={0,1,2,5,3,9,6,11,4,15,10,8,7,14,12,13};

