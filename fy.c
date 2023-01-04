#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

// sort's assorted
#include "global.h"
#include "fy.h"
typedef unsigned int d_type; // ソートするキーの型

static void merge1(d_type array[], d_type work[], int left, int mid, int right)
{
    d_type *array_i = array + left;
    d_type *array_j = array + mid;
    d_type *work_k = work + left;
    d_type *array_x = array + (mid - 1);
    d_type *array_y = array + (right - 1);
    d_type *work_z = work + (right - 1);
    int kosuu = right - left;

    for (int c = kosuu / 2; c > 0; c--) {
        MIN_TO_WORK  ;  MAX_TO_WORK
    }
    if (kosuu % 2) {
        //データ数が奇数のとき
        if (array_i == array_x)
            *work_k = *array_i;
        else
            *work_k = *array_j;
    }
}

void insertion_sortA(d_type array[G_N], int right)
{
    if (array[0] > array[1])               // 高速化のため
        SWAP(d_type, array[0], array[1]);

    for (int i = 2; i < right; i++) {
        int temp = array[i];
        if (array[i - 1] > temp) {
            int j = i;
            do {
                array[j] = array[j - 1];
                j--;
            } while (j > 0 && array[j - 1] > temp);
            array[j] = temp;
        }
    }
}

static void insertion_sortW(d_type *array, d_type *work, int right)
{
    if (array[0] <= array[1]) {   // 高速化のため
        work[0] = array[0];
        work[1] = array[1];
    } else {
        work[0] = array[1];
        work[1] = array[0];
    }

    for (int i = 2; i < right; i++) {
        if (work[i - 1] > array[i]) {
            int j = i;
            do {
                work[j] = work[j - 1];
                j--;
            } while (j > 0 && work[j - 1] > array[i]);
            work[j] = array[i];
        } else {
            work[i] = array[i];
        }
    }
}

/* merge_sortA() と merge_sortW() は相互再帰？ */
static void merge_sortA(d_type array[], d_type work[], int left, int right);
static void merge_sortW(d_type array[], d_type work[], int left, int right);

#define UTIKIRI 7
static void merge_sortA(d_type array[], d_type work[], int left, int right)
{
    int kosuu = right - left;

    if (kosuu <= UTIKIRI) {
        insertion_sortA(array + left, kosuu);
        return;
    }
    int mid = left + kosuu / 2;
    merge_sortW(array, work, left, mid);
    merge_sortW(array, work, mid, right);
    merge1(work, array, left, mid, right);
}

static void merge_sortW(d_type array[], d_type work[], int left, int right)
{
    int kosuu = right - left;

    if (kosuu <= UTIKIRI) {
        insertion_sortW(array + left, work + left, kosuu);
        return;
    }
    int mid = left + kosuu / 2;
    merge_sortA(array, work, left, mid);
    merge_sortA(array, work, mid, right);
    merge1(array, work, left, mid, right);
}

/* xをkビット右シフトし、その左iビットをとりだす */
static int bits(int x, int k, int i)
{
    return (x >> k) & ~(~0 << i);
}

typedef struct {
    int rand;           // 出席番号
    unsigned short ind; // 成績
} SRT;

/* 基数ソート  -- 再帰の仕方が quick sort() に似てる */
static void radix_sort(SRT array[], int left, int right, int bit)
{
    if (left < right && bit >= 0) {
        int i = left;
        int j = right;

        do {
            while (bits(array[i].rand, bit, 1) == 0 && (i < j))
                i++; // 左端から探索
            
            while (bits(array[j].rand, bit, 1) == 1 && (i < j))
                j--; // 右端から探索

            if (i != j) SWAP(SRT, array[i], array[j]);

        } while (i != j);

        if (bits(array[right].rand, bit, 1) == 0)
            j++;

        radix_sort(array, left, j - 1, bit - 1); // 0の部分をソート
        radix_sort(array, j, right, bit - 1);    // 1の部分をソート
    }
}

void merge_rand(unsigned short a[], int n)  /// n = G_N で呼ばれる、不要？
{
    SRT data[G_N];
    time_t t;

    srand(time(&t));

    // データ作成
    for (int i = 0; i < G_N; i++) {
        data[i].rand = random() % 0xffffffff;  /// 剰余不要？
        data[i].ind = i;
    }

    /* 直接基数法 */
    printf("straight_radix_sort\n");
    printf("Before sort: ");

    /* ソート後のデータを表示     */
    radix_sort(data, 0, G_N - 1, 32);

    for (int i = 0; i < G_N; i++) {
        //printf("%d,%d %d\n", i, data[i].ind, data[i].rand);
        a[i] = data[i].ind;
    }
}
