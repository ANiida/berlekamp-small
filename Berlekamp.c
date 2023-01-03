// date 20200331:xgcdの終了条件が２つになってしまう。ogcdとxgcdで使い分ける。
// date : 20200326 鍵生成が det と deta で高い確率で一致する。detaは並列処理。
// date 20200229 : pattarson algorithm implementation ver 1.0
//  xgcd & osqrtを追加した
// date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
// auther    : the queer who thinking about cryptographic future
// code name :  一変数多項式演算ライブラリのつもり
// status    : now in debugging (ver 0.8)
//  0ベクトルが出ないように生成多項式のトレースチェックを入れた。
// date      :  20160310,20210419
// auther    : the queer who thinking about cryptographic future
// code name : OVP - One Variable Polynomial library with OpenMP friendly
// status    : now in debugging (ver 0.9)
//
// 2022/06?? rubato 作業始める
// 2022/0701 rubato 大ナタ振るう

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>

#include <omp.h>
// #include "4096.h"
#include "global.h"
#include "struct.h"
#include "gf.h"
#include "val.h"

extern unsigned long xor128(void);
extern int mlt(int x, int y);
extern int mltn(int n, int a);

extern MTX inv_S;
extern MTX norm_S; /// extern MTX S;
extern MTX mulmat(MTX A, MTX B, int i);

void print_trace(void);
void random_shuffle(unsigned short *array, size_t size);
int mkRS(MTX cc, MTX *R);
int is_reg(MTX cc, MTX *R);
unsigned short merge_rand(unsigned short *a, int n); /// in fy.c

/* Goppa多項式 */
static unsigned short g[G_K + 1] = {0}; /// ginit(), ogt(), mkpol(), mkd()

//// static unsigned int AA = 0; 削除
//// static unsigned int B = 0;  削除
//// static MTX H = {0}; ==> pk_gen() 内に移動

// 有限体の元の逆数
static unsigned short oinv(unsigned short a)
{
    if (a == 0)
        return 0;

    for (unsigned short i = 0; i < G_N; i++)
    {
        if (gf[mlt(fg[a], i)] == 1)
            return i;
    }
    printf("no return \n");
    return 0; //  exit (1);
}

// aに何をかけたらbになるか
static unsigned short equ(unsigned short a, unsigned short b)
{
    int i;

    for (i = 0; i < G_N; i++)
    {
        if (gf[mlt(fg[a], fg[i])] == b)
            break;
    }
    return i;
}

// OP型からベクトル型への変換
static vec o2v(OP f)
{
    vec a = {0};

    for (int i = 0; i < G_DEG; i++)
    {
        if (f.t[i].a > 0 && f.t[i].n < G_DEG)
            a.x[f.t[i].n] = f.t[i].a;
    }
    return a;
}

// ベクトル型からOP型への変換
static OP v2o(vec a)
{
    OP f = {0};
    int j = 0;

    for (int i = 0; i < G_DEG; i++)
    {
        if (a.x[i] > 0)
        {
            f.t[j].n = i;
            f.t[j++].a = a.x[i];
        }
    }
    return f;
}

// 停止コマンド
#define LINESIZE 1 //// <= ヘン！　sizeof(line) ???
static void wait(void)
{
    char line[30];
    char *result;

    // 何か表示させたほうが良いだろう
    printf(" (enter number and hit return) ");
    fflush(stdout);

    //// ここで fgets() がエラーになることで、全体の動作がとまらないと思われる
    if ((result = fgets(line, LINESIZE, stdin)) != NULL)
        printf("The string is %s\n", result);
}

// OP型を正規化する
static OP conv(OP f)
{
    vec v = o2v(f);
    OP g = v2o(v);
    return g;
}

// 多項式の次数(default)
static int deg(vec a)
{
    int n = 0;
    bool flg = false;

    for (int i = 0; i < G_DEG; i++)
    {
        if (a.x[i] > 0)
        {
            n = i;
            flg = true;
        }
    }
    return flg ? n : 0; //// return n; ????
}

// 項の数
static int terms(OP f)
{
    int count = 0;
    for (int i = 0; i < G_DEG; i++)
        if (f.t[i].a > 0)
            count++;

    return count;
}

// 多項式の次数(degのOP型)
static int odeg(OP f)
{
    if (f.t[0].a == 0)
        return 0;

    int j = 0;
    for (int i = 0; i < G_DEG; i++)
    {
        if (j < f.t[i].n && f.t[i].a > 0)
            j = f.t[i].n;
    }
    return j;
}

// 多項式を表示する（OP型）
static void oprintpol(OP f)
{
    int n;

    f = conv(f);
    n = odeg(f);
    printf("n=%d\n", n);
    printf("terms=%d\n", terms(f));
    printf("deg=%d\n", odeg(f));

    for (int i = n; i > -1; i--)
    {
        if (f.t[i].a > 0)
            printf("%ux^%u+", f.t[i].a, f.t[i].n);
    }
}

static void op_print_raw(const OP f)
{
    puts("op_print_raw:");
    for (int i = 0; i < G_DEG; i++)
    {
        if (f.t[i].a > 0)
            printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
    }
}

static bool op_verify(const OP f)
{
    bool end = false;
    unsigned short n_max = 0;
    for (int i = 0; i < G_DEG; i++)
    {
        if (end && (f.t[i].n != 0 || f.t[i].a != 0))
        {
            op_print_raw(f);
            printf("found data after end: i=%d\n", i);
            print_trace();
            fflush(stdout);
            return false;
        }
        if (f.t[i].a == 0)
        {
            end = true;
            continue;
        }
        if (f.t[i].n + 1 <= n_max)
        {
            op_print_raw(f);
            printf("found invalid order: i=%d\n", i);
            print_trace();
            fflush(stdout);
            return false;
        }
        n_max = f.t[i].n + 1;
    }
    return true;
}

#if 0
#endif

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
static OP oadd(OP f, OP g)
{
    f = conv(f);
    g = conv(g);
    assert(op_verify(f));
    assert(op_verify(g));

    vec a = o2v(f);
    vec b = o2v(g);
    vec c = {0};

    for (int i = 0; i < G_DEG; i++)
    {
        c.x[i] = a.x[i] ^ b.x[i];
    }
    OP h = v2o(c);
    assert(op_verify(h));
    return h;
}

// 多項式を項ずつ掛ける
static OP oterml(OP f, oterm t)
{
    f = conv(f);
    assert(op_verify(f));

    OP h = {0};
    for (int i = 0; i < G_DEG; i++)
    {
        h.t[i].n = f.t[i].n + t.n;
        h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
    }

    h = conv(h);
    assert(op_verify(h));
    return h;
}

// 多項式の掛け算
static OP omul(OP f, OP g)
{
    f = conv(f);
    g = conv(g);
    assert(op_verify(f));
    assert(op_verify(g));

    int k = odeg(f);
    int l = odeg(g);
    if (k < l)
    {
        k = l;
    }

    OP h = {0};
    for (int i = 0; i <= k; i++)
    {
        oterm t = g.t[i];
        OP e = oterml(f, t);
        h = oadd(h, e);
    }
    assert(op_verify(h));
    return h;
}

// リーディグタームを抽出(default)
static oterm LT(OP f) //// leadingTerm() ???
{
    oterm t = {0};

    for (int i = 0; i < G_DEG; i++)
    {
        if (f.t[i].a > 0)
        {
            t.n = f.t[i].n;
            t.a = f.t[i].a;
        }
    }
    return t;
}

// 多項式を単行式で割る
static oterm LTdiv(OP f, oterm t)
{
    oterm tt, s = {0};

    tt = LT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = gf[mlt(fg[tt.a], oinv(t.a))];
        s.n = tt.n;
    }

    return s;
}

// 多項式を表示する(default)
static void printpol(vec a)
{
    int n = deg(a);
    assert(n >= 0);
    for (int i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            // if(a.x[i]>1)
            printf("%u*", a.x[i]);
            if (i > 0)
                printf("x^%d", i);
            // if (i > 0)
            if (i > 0)
                printf("+");
        }
    }
}

// 多項式の剰余を取る
static OP omod(OP f, OP g)
{
    int n = LT(g).n;

    if (LT(f).n < n)
        return f;

    oterm b = LT(g);
    assert(b.a != 0 && b.n != 0);
    while (LT(f).n > 0 && LT(g).n > 0)
    {
        oterm c = LTdiv(f, b);
        OP h = oterml(g, c);
        f = oadd(f, h);
        if (odeg((f)) == 0 || odeg((g)) == 0)
            break;

        if (c.n == 0 || b.n == 0)
            break;
    }
    return f;
}

// 多項式のべき乗余
static OP opowmod(OP f, OP mod, int n)
{
    // printpol(o2v(mod));
    // printf(" =mod %d\n",LT(mod).n);
    // 繰り返し２乗法
    for (int i = 1; i < n + 1; i++)
        f = omod(omul(f, f), mod);

    return f;
}

// 多項式の代入値
static unsigned short trace(OP f, unsigned short x)
{
    unsigned short u = 0;
    int d = deg(o2v(f));

    for (int i = 0; i < d + 1; i++)
    {
        u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
    }
    return u;
}

static unsigned short v2a(oterm a)
{
    if (a.a == 0)
        return 0;

    for (int j = 0; j < G_M; j++)
    {
        if (gf[j] == a.a && a.a > 0)
        {
            return j - 1;
        }
    }
    return 0;
}

static void printsage(vec a)
{
    printf("poly=");
    for (int i = 0; i < G_DEG; i++)
    {
        if (a.x[i] > 0)
        {
            oterm b;

            b.a = a.x[i];
            b.n = i;

            int j = v2a(b);
            printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
        }
    }
}

static OP gcd(OP a, OP b)
{
    OP r = {0}, h = {0}; // , tmp = {0};

    h.t[0].a = 1;
    h.t[0].n = 0;

    if (odeg(a) < odeg(b))
    {
        OP tmp = a;
        a = b;
        b = tmp;
    }
    /* 自然数 a > b を確認・入替 */
    if (odeg(a) < odeg(b))
    {
        OP tmp = a;
        a = b;
        b = tmp;
    }

    r = omod(a, b);
    while (odeg(r) > 0)
    {
        a = b;
        b = r;
        r = omod(a, b);
        if (LT(r).a == 0)
            return b;
    }

    return (LT(r).a == 0) ? b : h;
}

#if 0
#endif
#if 0
#endif
#if 0
#endif
#if 0
#endif
#if 0
#endif
#if 0
#endif

static OP kof(unsigned short c, OP f)
{
    int i, k;
    vec b = {0}, h = {0};
    OP g = {0};

    c = fg[c];
    b = o2v(f);
    k = deg(b);
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
    }
    g = v2o(h);

    return g;
}

// ランダム多項式の生成
static void ginit(void)
{
    int j, count = 0, k;
    unsigned short gg[G_K + 1] = {0};

    printf("in ginit\n");

    memset(g, 0, sizeof(g)); /// mkpol() から移動
    g[G_K] = 1;              // xor128();
    g[0] = 1;                // random() % G_N; // or G_N

    k = random() % (G_K - 1);
    if (k > 0)
    {
        while (count < k)
        {
            printf("in whule\n");
            j = random() % (G_K);
            if (j < G_K && j > 0 && g[j] == 0)
            {
                g[j] = random() % G_M;
                count++;
            }
        }
    }

    for (j = 0; j < G_K + 1; j++)
        gg[j] = g[G_K - j]; /// 逆順にする？

    memcpy(g, gg, sizeof(g));
}

// 整数からベクトル型への変換
static vec i2v(unsigned int n)
{
    vec v = {0};
    int i = 0;

    while (n > 0)
    {
        v.x[i++] = n % 2;
        n = (n >> 1);
    }
    return v;
}

// ベクトル型から整数への変換
static unsigned int v2i(vec v)
{
    unsigned int d = 0, e = 0;

    for (int i = 0; i < deg(v) + 1; i++)
    {
        e = v.x[i];
        d ^= (e << i);
    }
    return d;
}

// 配列からベクトル表現の多項式へ変換する
static vec Setvec(int n)
{
    vec v = {0};

    for (int i = 0; i < n; i++)
        v.x[n - 1 - i] = c[i];

    return v;
}

// chen探索
static vec chen(OP f)
{
    vec e = {0};
    int n = odeg(f);
    int count = 0;

    for (int x = 0; x < G_N; x++)
    {
        unsigned short z = 0;

        for (int i = 0; i <= n; i++)
        {
            if (f.t[i].a > 0)
                z ^= gf[mlt(mltn(f.t[i].n, fg[x]), fg[f.t[i].a])];
        }
        if (z == 0)
        {
            e.x[count] = x;
            count++;
            printf("%d\n", x);
        }
    }

    return e;
}
/*
// GF(2^m) then set m in this function.
static int ben_or(OP f)
{
    int i, n, flg = 0;
    OP s = {0}, u = {0}, r = {0};
    vec v = {0};
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    int m = G_E;
    // m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = v2o(v);
    r = s;
    n = deg(o2v(f));
    printf("n=%d\n", n);

    if (n == 0)
        return -1;

    // r(x)^{q^i} square pow mod
    i = 0;
    while (i < n / 2) {
        printf("iii=%d\n", i);
        flg = 1;
        // irreducible over GH(8192) 2^13
        r = opowmod(r, f, m);

        // irreducible over GF2  // r=omod(opow(r,2),f);
        u = oadd(r, s);
        if (deg(o2v(u)) == 0 && LT(u).a == 0)
            return -1;

        if (deg(o2v(u)) == 0 && LT(u).a == 1) {
            i++;
            flg = 0;
        }
        if (deg(o2v(u)) > 0)
            u = gcd(f, u);

        if (deg(o2v(u)) > 0)
            return -1;

        if (flg == 1)
            i++;
    }

    return 0;
}
*/
// 配列の値を係数として多項式に設定する
static OP setpol(unsigned short f[], int n)
{
    OP g;
    vec a;

    memset(c, 0, sizeof(c));
    memcpy(c, f, 2 * n);
    a = Setvec(n);
    g = v2o(a);
    return g;
}

// バイナリ型パリティチェック行列を生成する
static MTX bdet()
{
    MTX R = {0};

    for (int i = 0; i < G_N; i++)
    {
        for (int j = 0; j < G_K; j++)
        {
            int l = mat[i][j];

            for (int k = 0; k < G_E; k++)
            {
                R.x[i][j * G_E + k] = l % 2;
                l = (l >> 1);
            }
        }
    }

    for (int i = 0; i < G_N; i++)
    {
        for (int j = 0; j < G_DEG; j++)
        {
            printf("%d,", R.x[i][j]);
        }
        printf("\n");
    }
    return R;
}

static MTX bd2()
{
    vec v = {0};
    MTX R = {0};

    for (int i = 0; i < G_N; i++)
    {
        for (int j = 0; j < G_K21; j++)
        {
            int l = bm[i][j];
            printf("bm==%d %d\n", l, j);

            v = i2v(l);

            for (int k = 0; k < G_E; k++)
            {
                R.x[i][j * G_E + k] = v.x[k];
            }
        }
    }

    return R;
}

static unsigned short HH[G_N][G_K];

static MTX toByte(MTX SH, int kk)
{
    vec v = {0};
    MTX R = {0};
    int i, j, k, cnt;

    memset(HH, 0, sizeof(HH));
    printf("HH=");

    for (i = 0; i < G_N; i++)
    {
        printf("%d\n", i);
        for (j = 0; j < kk; j++)
        {
            cnt = 0;
            for (k = j * G_E; k < j * G_E + G_E; k++)
                v.x[cnt++] = SH.x[i][k];

            R.x[i][j] = v2i(v);
        }
    }
    printf("end of byte\n");
    return R;
}

static vec b2v(vec v)
{
    vec x = {0}, z = {0};

    for (int j = 0; j < G_K; j++)
    {
        int i = 0;
        for (int k = j * G_E; k < j * G_E + G_E; k++)
            x.x[i++] = v.x[k];
        z.x[j] = v2i(x);
    }
    for (int i = 0; i < G_K; i++)
        printf("%d,", z.x[i]);

    return z;
}

// 秘密置換を生成する
static void Pgen()
{
    CNT++;
    memset(P, 0, sizeof(P));
    merge_rand(P, G_N);

    for (int i = 0; i < G_N; i++)
        inv_P[P[i]] = i;
}

OP synd(unsigned short zz[], int kk)
{
    unsigned short syn[G_K] = {0}, s = 0;
    OP f = {0};

    printf("in synd2\n");

    for (int i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        for (int j = 0; j < G_M; j++)
        {
            s ^= gf[mlt(fg[zz[j]], fg[mat[j][i]])];
        }
        syn[i] = s;
    }

    f = setpol(syn, kk);
    printpol(o2v(f));
    printf(" syn=============\n");
    return f;
}

static vec sina(unsigned short zz[], MTX R)
{
    vec v = {0}, z = {0};

    for (int i = 0; i < G_N; i++)
    {
        if (zz[i] > 0)
        {
            for (int j = 0; j < G_DEG; j++)
            {
                v.x[j] ^= R.x[i][j];
                printf("%d,", R.x[i][j]);
            }
        }
        printf("\n");
    }

    z = b2v(v);
    return z;
}

static unsigned short vb[G_K * 2][G_N] = {0};
static unsigned short gt[G_K * 2][G_K * 2] = {0};

static void van(int kk)
{
    printf("van der\n");

    for (int i = 0; i < G_N; i++)
        vb[0][i] = 1;

    for (int i = 1; i < kk; i++)
    {
        for (int j = 0; j < G_N; j++)
        {
            vb[i][j] = gf[mltn(i, fg[j])];
            printf("%d,", vb[i][j]);
        }
        printf("\n");
    }
}

void ogt(unsigned short pp[], int kk)
{
    int i, j;                          //// ←消すな、このまま
#pragma omp parallel for private(i, j) //// ←消すな

    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk - i; j++)
        {
            gt[i][j + i] = g[j];
        }
    }
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk; j++)
            printf("%d,", gt[i][j]);
        printf("\n");
    }
}

static OP mkpol() //// mkpol2() との違いは？？
{
    int j, ii = 0;
    OP w = {0};

    do
    {
        int k, flg;

        j = k = flg = 0;
        //  memset(g, 0, sizeof(g));  ==> ginit() へ移動
        memset(w.t, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > 100)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (int i = 0; i < G_K; i++)
        {
            if (g[G_K - 1] > 0)
                flg = 1;
            if (i % 2 == 1 && g[i] > 0 && i < G_K)
                k++;
        }

        // 偶数項だけにならないようにする
        if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        {
            w = setpol(g, G_K + 1);
            j = 1;
        }

    } while (j == 0);

    printpol(o2v(w));
    printf(" ==g\n");

    return w;
}

static OP mkd(OP w, int kk)
{
    int i, j, ii, ll = -1;
    unsigned short ta[G_N] = {0};
    vec v = {0};
    OP r = {0};

aa:
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。
    //  デフォルトでGF(8192)    //既約多項式しか使わない。

    ii = 0;
    // irreducible goppa code (既役多項式が必要なら、
    // ここのコメントを外すこと。)  ????

    vec pp = {0};
    vec tt = {0};
    int l = -1;
    while (l < 0)
    {
        for (i = 0; i < G_K; i++)
            pp.x[i] = rand() % G_N;
        mykey(tt.x, pp);
        tt.x[G_K] = 1;
        l = ben_or(tt);
        if (l == 0)
        {
            printf("\n");
            printsage(tt);
            printf(" ==irr\n");
            break;
        }
        if (l == 0)
            break;
    }
    // exit(1);

    memset(ta, 0, sizeof(ta));
    printpol((tt));
    r = v2o(tt);
    // 多項式の値が0でないことを確認
    for (i = 0; i < G_N; i++)
    {
        ta[i] = trace(r, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            goto aa;
        }
    }

    // 多項式を固定したい場合コメントアウトする。
    oprintpol(r);
    printf("\n");
    printsage(o2v(r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    van(kk);

#if 0
    /// 必要か？
    memset(g, 0, sizeof(g));  g[0] = 1;
#endif
    ogt(g, kk);

    printf("\nすげ、オレもうイキそ・・・\n");
    for (i = 0; i < G_N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = vb[j][i];
        }
    }

    return w;
}

static void half(int kk)
{
    for (int i = 0; i < G_N; i++)
        bm[i][0] = mat[i][0];

    for (int i = 1; i < kk; i++)
    {
        for (int j = 0; j < G_N; j++)
            bm[j][i] = mat[j][i * 2 - 1];
    }
}

// Niederreiter暗号の公開鍵を作る(RS)
static MTX mk_pub()
{
    OP w = mkd(w, G_K * 2);
    MTX g_bin = {0};

    half(G_K21);

    oprintpol(w);
    printf("\n");
    printsage(o2v(w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    MTX r = bd2();
    printf("deco_rev= ");
    for (int i = 0; i < G_AY; i++)
        printf("%d,", r.x[1][i] ^ r.x[2][i] ^ r.x[3][i]);
    printf("\n");
    Pgen();

    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(norm_S.x, 0, sizeof(norm_S.x));
        for (int i = 0; i < G_AY; i++)
        {
            for (int j = 0; j < G_AY; j++)
                norm_S.x[i][j] = xor128() % 2;
        }
    } while (mkRS(norm_S, &inv_S) == -1);

    MTX z = mulmat(norm_S, r, 2);
    printf("Zz=\n");
    for (int j = 0; j < G_N; j++)
    {
        for (int i = 0; i < G_AY; i++)
        {
            g_bin.x[j][i] = z.x[P[j]][i];
        }
    }

    return toByte(g_bin, G_K21);
}

// Niederreiter暗号の公開鍵を作る(Goppa)
static MTX pk_gen()
{
    OP w = mkd(w, G_K);

    oprintpol(w);
    printf("\n");
    printsage(o2v(w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    Pgen();
    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(norm_S.x, 0, sizeof(norm_S.x));
        for (int i = 0; i < G_DEG; i++)
        {
            for (int j = 0; j < G_DEG; j++)
                norm_S.x[i][j] = xor128() % 2;
        }
    } while (is_reg(norm_S, &inv_S) == -1);

    MTX Q = bdet();
    MTX H = mulmat(norm_S, Q, 1);
    MTX O_bin = {0};
    for (int i = 0; i < G_DEG; i++)
    {
        for (int j = 0; j < G_N; j++)
        {
            O_bin.x[j][i] = H.x[P[j]][i];
        }
    }
    return O_bin;
}

static OP dec(unsigned short ss[])
{
    unsigned int ch[G_DEG] = {0};
    unsigned char h2o[G_DEG] = {0};
    vec v;

    printf("!1\n");
    for (int i = 0; i < G_K; i++)
    {
        v = i2v(ss[i]);
        for (int j = 0; j < G_E; j++)
            ch[i * G_E + j] = v.x[j];
    }
    for (int i = 0; i < G_DEG; i++)
        printf("%d", ch[i]);
    printf("\n");

    for (int i = 0; i < G_DEG; i++)
    {
        for (int j = 0; j < G_DEG; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }

    unsigned short uk[G_K];
    for (int i = 0; i < G_K; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        for (int j = 0; j < G_E; j++)
            v.x[j] = h2o[i * G_E + j];
        uk[i] = v2i(v);
    }
    for (int i = 0; i < G_K; i++)
        printf("%d,", uk[i]);
    printf("\n");

    return setpol(uk, G_K);
}

static int ero2(vec v)
{
    unsigned short ya[G_N] = {0}, xa[G_N] = {0};
    int count = 0;

    for (int i = 0; i < G_T; i++)
    {
        if (i == 0)
        {
            xa[v.x[i]] = 1;
            count++;
        }
        if (i > 0 && v.x[i] > 0)
        {
            xa[v.x[i]] = 1;
            count++;
        }
        if (i > 0 && v.x[i] == 0)
        {
            printf("baka %d %d\n", i, v.x[i]);
            printf("v.x[K-1]=%d\n", v.x[G_K - 1]);
            break;
        }
    }

    int cnt = 0;
    for (int i = 0; i < G_N; i++)
        ya[i] = xa[P[i]];

    for (int i = 0; i < G_N; i++)
    {
        if (ya[i] > 0 && i == 0)
        {
            printf("error position=%d う\n", i);
            cnt = 1;
        }
        else if (ya[i] > 0)
        {
            if (cnt == 0)
            {
                printf("error position=%d う\n", i);
                cnt = 1;
            }
            else
            {
                printf("error position=%d お\n", i);
            }
        }
    }

    if (count == G_T)
    {
        printf("err=%dっ!! \n", count);
        /// B++; 削除
    }
    if (count < G_T)
    {
        printf("error is too few\n");
        /// AA++; 削除
        printf("へげえええーっ\n");
        exit(1);
    }

    return count;
}

static void mkerr(unsigned short *z1, int num)
{
    int j, l;

    j = 0;

    memset(z1, 0, sizeof(*z1));

    while (j < num)
    {
        l = random() % G_N;
        // printf ("l=%d\n", l);
        if (0 == z1[l])
        {
            z1[l] = 1;
            // printf("l=%d\n", l);
            j++;
        }
    }
}

static vec newhalf(unsigned short e[G_K21])
{
    int i, k;
    vec v = {0};
    unsigned short t[G_K] = {0};

    for (i = 0; i < G_K21; i++)
    {
        t[i] = e[i];
        printf("e=%d\n", e[i]);
    }

    for (i = 0; i < G_K21; i++)
    {
        printf("t=%d\n", t[i]);
    }

    v.x[0] = t[0];
    v.x[1] = t[1];
    k = 2;
    for (i = 2; i < G_K; i++)
    {
        if (i % 2 == 1)
        {
            v.x[i] = t[k];
            k++;
        }

        if (i % 2 == 0)
            v.x[i] = gf[mlt(fg[v.x[i / 2]], fg[v.x[i / 2]])];
    }

    return v;
}

static vec bfd(unsigned short ss[])
{
    unsigned int ch[G_DEG * 2] = {0};
    int count = 0;
    vec v;

    for (int i = 0; i < G_K21; i++)
    { // count=(K/2+1)*E-1;
        v = i2v(ss[i]);
        for (int j = 0; j < G_E; j++)
        {
            ch[count] = v.x[j];
            count++;
        }
    }
    printf("bfd_bin=\n");
    for (int i = 0; i < G_AY; i++)
        printf("%d", ch[i]);
    printf("\n");

    unsigned char h2o[G_DEG * 2] = {0};
    for (int i = 0; i < G_AY; i++)
    {
        for (int j = 0; j < G_AY; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }

    printf("deco_bin=\n");
    for (int i = 0; i < G_AY; i++)
        printf("%d,", h2o[i]);
    printf("\n");

    unsigned short uk[G_K] = {0};
    count = 0;
    for (int i = 0; i < G_K21; i++)
    { // count=(K/2+1)*E-1;
        memset(v.x, 0, sizeof(v.x));
        for (int j = 0; j < G_E; j++)
        {
            v.x[j] = h2o[count];
            count++;
        }
        uk[i] = v2i(v);
    }

    printf("bm_int=\n");
    for (int i = 0; i < G_K21; i++)
        printf("%d,", bm[1][i] ^ bm[2][i] ^ bm[3][i]);
    printf("\n");

    printf("u_int=\n");
    for (int i = 0; i < G_K21; i++)
        printf("%d,", uk[i]);
    printf("\n");

    return o2v(setpol(uk, G_K21));
}

/*
 * Berlekamp-Massey Algorithm
 */
static OP bma(unsigned short s[], int kk)
{
    int j, k, ll, l;
    int d[2 * G_K + 1] = {0};
    OP lo[2 * G_K + 1] = {0};
    OP b[2 * G_K + 1] = {0};
    OP t[2 * G_K + 1] = {0};
    OP h = {0};
    OP g = {0};
    vec v = {0};
    vec x = {0};

    x.x[1] = 1;
    h = v2o(x);
    v.x[0] = 1;

    lo[0] = v2o(v);
    b[0] = lo[0];
    ll = 0;
    for (j = 1; j < G_T * 2 + 1; j++)
    {
        v = o2v(lo[j - 1]);
        k = 0;

        l = deg(o2v(lo[j - 1]));
        for (int i = 1; i < l + 1; i++)
        {
            k ^= gf[mlt(fg[v.x[i]], fg[s[j - i]])];
        }
        d[j] = s[j] ^ k;

        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = omul(b[j - 1], h);
        }
        else
        {
            g = omul(kof(d[j], h), b[j - 1]);
            t[j] = oadd(lo[j - 1], g);
            lo[j] = t[j];
            if (ll * 2 > (j - 1))
            {
                b[j] = omul(b[j - 1], h);
            }
            else
            {
                b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
                ll = j - ll;

                if (j == 2 * G_T)
                {
                    if (!(d[G_T * 2 - 1] == 0 &&
                          d[G_T * 2 - 3] == 0 &&
                          odeg(lo[j - 1]) == G_T) ||
                        !(odeg(lo[j - 1]) == G_T))
                    {
                        if ((d[G_T * 2 - 1] == 0 && odeg(lo[j - 2]) == G_T - 1))
                        {
                            lo[j - 1] = omul(lo[j - 2], h);
                        }
                    }
                    break;
                }
            }
        }
        printf("l=%d\n", ll);
        k = 0;
    }

    k = 0;
    if (odeg(lo[j - 1]) == G_T)
    {
        x = chen(lo[j - 1]);
    }
    else
    {
        printf("baka==\n");
        exit(1);
    }

    int count = 0;
    for (int i = 0; i < deg(x) + 1; i++)
    {
        if (x.x[i] >= 0)
        {
            printf("xx[%d]=1\n", (x.x[i]));
            count++;
        }

        if (x.x[i] == 0)
            k++;

        if (k > 1)
        {
            printf("baka0\n");
            exit(1);
        }
    }

    if (count < G_T)
    {
        printf("vaka in bms %d\n", count);
        exit(1);
    }

    return lo[j - 1]; // return count;
}

static OP sendrier2(unsigned short zz[G_N], MTX L)
{
    unsigned short s[G_K + 1] = {0}, rt[G_K21] = {0};
    int j;

    for (j = 0; j < G_N; j++)
    {
        if (zz[j] > 0)
        {
            for (int k = 0; k < G_K21; k++)
            {
                rt[k] ^= L.x[j][k];
            }
        }
    }
    for (int i = 0; i < G_K21; i++)
        printf("%d,", rt[i]);
    printf("\n");

    vec u = {0};
    vec x = bfd(rt);
    for (int i = 0; i < G_K21; i++)
    {
        u.x[G_K / 2 - i] = x.x[i];
    }
    for (int i = 0; i < G_K21; i++)
        printf("%d,%d\n", u.x[i], x.x[i]);
    printf("\n");

    printf("rt=\n");
    for (int i = 0; i < G_K21; i++)
        printf("%d,", u.x[i]);
    printf("\n");
    printf("bm_in se2 == %d || ", j);

    vec v = newhalf(u.x);

    printf("P= ");
    for (int i = 0; i < G_N; i++)
        printf("%d,", P[i]);
    printf("\n");

    memset(s, 0, sizeof(s));
    for (int i = 0; i < G_K; i++)
        s[i + 1] = v.x[i];

    printf("rt_deco= ");
    return bma(s, G_K);
}

// 言わずもがな
int main(void)
{
    int i;
    unsigned short zz[G_N] = {0};
    OP f = {0}, r = {0};
    vec v = {0}, x = {0};
    MTX R = {0}, O = {0};
    unsigned short s[G_K + 1] = {0};
    struct timespec ts;

    if (timespec_get(&ts, TIME_UTC) == 0)
    {
        exit(1);
    }
    srandom(ts.tv_nsec ^ ts.tv_sec); /* PRNG にシードを設定する */
    if (G_K > G_N)
        printf("configuration error! G_K is bigger than G_N\n");

    memset(mat, 0, sizeof(mat));

    // 公開鍵を生成する(Niederreiterとは異なる)
    // 鍵サイズ G_K : Goppa Code
    R = pk_gen();

    // エラーベクトルの初期化
    memset(zz, 0, sizeof(zz));
    // 重み G_T のエラーベクトルを生成する
    mkerr(zz, G_T);
    // 暗号文の生成(s=eH)
    // x = sin2(zz, R);
    x = sina(zz, R);
    // 復号化１(m'=sS^{-1})
    r = dec(x.x);
    v = o2v(r);
    for (i = 0; i < G_K; i++)
        s[i + 1] = v.x[i];

    // Berlekamp-Massey Algorithm
    f = bma(s, G_K);
    x = chen(f);
    // 平文の表示(m=m'P^{-1})
    ero2(x);
    for (i = 0; i < G_N; i++)
        if (zz[i] > 0)
            printf("err=%d\n", i);
    wait();
    exit(1);

    // debugging
    O = mk_pub(); // 鍵サイズ(K/2 Reed-Solomon)
    memset(zz, 0, sizeof(zz));
    for (i = 0; i < G_T; i++)
        zz[i] = 1;

    f = sendrier2(zz, O);
    x = chen(f);
    ero2(x);
    printf("aaa\n");

    for (i = 0; i < G_N; i++)
        if (zz[i] > 0)
            printf("i= %d,", i);
    printf("\n");

    if (odeg(f) < G_T)
    {
        printpol(o2v(r));
        printf("==goppa\n");
        for (i = 0; i < G_N; i++)
            printf("%d,", zz[i]);
        printf("\n");
        exit(1);
    }

    return 0;
}
