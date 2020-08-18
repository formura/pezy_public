#pragma once

// 初期条件の選択
// vortex問題をやるときは VORTEX をdefineすること
/* #define VORTEX */

// NX,NY,NZ は床のサイズ
// NX = NY = NZ でなければならない
#define NX 44
#define NY 44
#define NZ 44

// NTはTemporal Blockingの段数
#define NT 16

// MX,MY,MZ は床の数
#define MX 20
#define MY 20
#define MZ 20

// Nsはスキームの袖サイズ
#define Ns 2

#define T_MAX 50
#define FILTER_INTERVAL 10

#define LX (NX * MX - 2 * Ns * NT)
#define LY (NY * MY - 2 * Ns * NT)
#define LZ (NZ * MZ - 2 * Ns * NT)

// 乱流のエネルギースペクトルのパラメータ
// E(k) = A*k^4*exp(-B*k^2)
#define A 0.25
#define B 0.54

#define FLOP_PER_ELEM 2829

constexpr double CONST_h = 1.0 / LX;
constexpr double CONST_cfl = 1.0;
constexpr double CONST_dt = CONST_cfl * CONST_h * CONST_h;
constexpr double CONST_re = 150.0;
constexpr double CONST_mu = 1.0;
constexpr double CONST_gm = 1.4;
constexpr double CONST_c = CONST_mu / CONST_re;
constexpr double CONST_c2 = (CONST_gm - 1) * CONST_mu / CONST_re;

#if !defined(__PZC_KERNEL__)
constexpr double h = CONST_h;
constexpr double cfl = CONST_cfl;
constexpr double dt = CONST_dt;
constexpr double re = CONST_re;
constexpr double mu = CONST_mu;
constexpr double gm = CONST_gm;
constexpr double c = CONST_c;
constexpr double c2 = CONST_c2;
#else

struct constants {
    double h = CONST_h;
    double cfl = CONST_cfl;
    double dt = CONST_dt;
    double re = CONST_re;
    double mu = CONST_mu;
    double gm = CONST_gm;
    double c = CONST_c;
    double c2 = CONST_c2;

    double num1 = 1.0;
    double num2 = 2.0;
    double num3 = 3.0;
    double num4 = 4.0;
    double num5 = 5.0;
    double num6 = 6.0;
    double num8 = 8.0;
    double num12 = 12.0;
    double num16 = 16.0;
    double num30 = 30.0;
    double num64 = 64.0;
    double num128 = 128.0;
    double num240 = 240.0;
    double num256 = 256.0;
    double num480 = 488.0;
    double num512 = 512.0;
    double num900 = 900.0;
    double num1024 = 1024.0;
    double num1920 = 1920.0;

    double rev1h4 = 1.0 / (1.0 * CONST_h * CONST_h * CONST_h * CONST_h);
    double rev2h3 = 1.0 / (2.0 * CONST_h * CONST_h * CONST_h);
    double rev12h1 = 1.0 / (12.0 * CONST_h);
    double rev12h2 = 1.0 / (12.0 * CONST_h * CONST_h);
    double rev24h4 = 1.0 / (24.0 * CONST_h * CONST_h * CONST_h * CONST_h);
    double rev144h2 = 1.0 / (144.0 * CONST_h * CONST_h);
    double rev144h3 = 1.0 / (144.0 * CONST_h * CONST_h * CONST_h);
    double rev144h4 = 1.0 / (2.0 * CONST_h * CONST_h * CONST_h * CONST_h);
    double rev1728h3 = 1.0 / (1728.0 * CONST_h * CONST_h * CONST_h);
    double rev1728h4 = 1.0 / (1728.0 * CONST_h * CONST_h * CONST_h * CONST_h);

    double num_1_3 = 1.0 / 3.0;
    double num_4_3 = 4.0 / 3.0;
};

#endif

template <typename T>
struct elem {
    T r, u, v, w, p;

    elem& operator=(const elem& r)
    {
        this->r = r.r;
        this->u = r.u;
        this->v = r.v;
        this->w = r.w;
        this->p = r.p;
        return *this;
    }

    // template <typename U>
    // elem& operator=(const elem<U>& r)
    // {
    //     this->r = r.r;
    //     this->u = r.u;
    //     this->v = r.v;
    //     this->w = r.w;
    //     this->p = r.p;
    //     return *this;
    // }
};

template <typename T>
__attribute__((always_inline)) elem<T> operator+(const elem<T>& l, const elem<T>& r)
{
    return elem<T>{ l.r + r.r, l.u + r.u, l.v + r.v, l.w + r.w, l.p + r.p };
}

template <typename T>
__attribute__((always_inline)) elem<T> operator-(const elem<T>& l, const elem<T>& r)
{
    return elem<T>{ l.r - r.r, l.u - r.u, l.v - r.v, l.w - r.w, l.p - r.p };
}

template <typename T, typename U>
__attribute__((always_inline)) elem<T> operator*(U d, const elem<T>& r)
{
    return elem<T>{ d * r.r, d * r.u, d * r.v, d * r.w, d * r.p };
}

template <typename T, typename U>
__attribute__((always_inline)) elem<T> operator*(const elem<T>& r, U d)
{
    return elem<T>{ d * r.r, d * r.u, d * r.v, d * r.w, d * r.p };
}

template <typename T>
__attribute__((always_inline)) elem<T> operator-(const elem<T>& r)
{
    return elem<T>{ -r.r, -r.u, -r.v, -r.w, -r.p };
}

template <typename T>
__attribute__((always_inline)) elem<T> operator+(const elem<T>& r)
{
    return r;
}

template <typename T>
struct elem2 {
    elem<T> p, h;
};

template <typename T>
struct walls {
    elem2<T> x[NT][NZ + 2 * Ns][NY + 2 * Ns][2 * Ns];
    elem2<T> y[NT][NZ + 2 * Ns][2 * Ns][NX + 2 * Ns];
    elem2<T> z[NT][2 * Ns][NY + 2 * Ns][NX + 2 * Ns];
};

template <typename T>
struct tmp_walls {
    elem2<T> x[MY][MZ][NT][NZ + 2 * Ns][NY + 2 * Ns][2 * Ns];
    elem2<T> y[MZ][MX][NT][NZ + 2 * Ns][2 * Ns][NX + 2 * Ns];
    elem2<T> z[MX][MY][NT][2 * Ns][NY + 2 * Ns][NX + 2 * Ns];
};

template <typename T>
struct floors {
    elem2<T> v[MZ * NZ][MY * NY][MX * NX];
};

template <typename T>
struct local_buf {
    elem<T> p[NZ + 2 * Ns][NY + 2 * Ns][NX + 2 * Ns];
    elem<T> h[NZ + 2 * Ns][NY + 2 * Ns][NX + 2 * Ns];
};

typedef struct
{
    int time_step;
    int lower_x;
    int lower_y;
    int lower_z;
    int upper_x;
    int upper_y;
    int upper_z;
    int offset_x;
    int offset_y;
    int offset_z;
} navi;

// void lowpass_filter(navi &, state &);

// void init(navi &, state &);
// void next(navi &, state &);
