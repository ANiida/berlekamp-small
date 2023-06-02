#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "global.h"

#define str_length 128
#define password_length 256

char password[password_length + 1];

unsigned long xor128(void)
{
  unsigned int a = 0;

  static unsigned long x = 123456789, y = 362436069, z = 521288629, w = 88675123;
  unsigned long t;

  a = rand();
  t = x ^ (a << 11);
  a = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

void seed(void)
{

  // * 変数宣言

  char str[str_length + 1];
  time_t t;
  int i, j, k, rnd;

  // * 乱数の初期化

  srand(clock() + time(&t));

  // * 乱数生成とパスワードの生成

  for (i = 0; i < str_length; i++)
  {
    for (j = 0; j < 2; j++)
    {
      k = i * 2 + j;
      do
      {
        rnd = rand();
        password[k] = str[i] + rnd;
      } while (!isalnum(password[k]));
    }
  }

  // NULL文字の挿入

  password[password_length] = '\0';

  // * パスワードの出力

  //    printf("生成パスワード：%s",password);

  return;
}

int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (N - 1)) + 1;
}

int mltn2(int n, int x)
{
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  return mlt(x, mltn2(n - 1, x));
}

int mltn(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = mlt(i, x);

  return i;
}

void rondom_permtation(unsigned char *a)
{
  int i, j, x;

  //  srand(clock() + time(&t));

  for (i = 0; i < K; i++)
  {
    a[i] = i;
  }
  for (i = 0; i < K - 2; i++)
  {
    // rand from i+1 to F-1
    j = (rand() % (K - 1 - i)) + i + 1;

    // swap a[i] and a[j]
    x = a[j];
    a[j] = a[i];
    a[i] = x;
  }
  if (a[K - 1] == K - 1)
  {
    a[K - 1] = a[K - 2];
    a[K - 2] = K - 1;
  }
}

int print_uint128(__uint128_t n)
{
  if (n == 0)
    return printf("0\n");

  char str[40] = {0};              // log10(1 << 128) + '\0'
  char *s = str + sizeof(str) - 1; // start at the end
  while (n != 0)
  {
    if (s == str)
      return -1; // never happens

    *--s = "0123456789"[n % 10]; // save last digit
    n /= 10;                     // drop it
  }
  return printf("%s", s);
}
