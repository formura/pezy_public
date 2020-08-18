// TODO: int to double変換命令が出ないようにする
// TODO: 変数受け渡しちゃんとする
// TODO: 前の実装では壁バッファ使ってないから消せるはずだから消す

#ifndef __PZC_KERNEL__
#define __PZC_KERNEL__
#endif

#include "config.h"
template <typename T>
void swap(T& a, T& b)
{
    T t = a;
    a = b;
    b = t;
}

template <typename _Tp, size_t _Nm>
struct array {
    _Tp data[3];

    _Tp& operator[](size_t __n) noexcept { return data[__n]; }
    const _Tp& operator[](size_t __n) const noexcept { return data[__n]; }
};

template <typename T, size_t X, size_t Y, size_t Z>
INLINE T& get_elem(T (&arr)[Z][X][Y], size_t x, size_t y, size_t z)
{
    return arr[z][x][y];
}

// FIXME: refactoring

template <typename T,
    size_t X1,
    size_t Y1,
    size_t Z1,
    size_t X2,
    size_t Y2,
    size_t Z2>
INLINE void copy_cubo(T (&dst)[X1][Y1][Z1],
    const T (&src)[X2][Y2][Z2],
    const array<int, 3>& dst_offset,
    const array<int, 3>& src_offset,
    const array<int, 3>& region)
{
    int gid = get_pid() * THREADS_PER_PE + get_tid();
    int gws = get_maxpid() * THREADS_PER_PE;

    for (int i = gid; i < region[0] * region[1] * region[2]; i += gws) {
        int x = i / region[2] / region[1];
        int y = i / region[2] % region[1];
        int z = i % region[2];

        get_elem(dst, x + dst_offset[0], y + dst_offset[1], z + dst_offset[2])
            = get_elem(src, x + src_offset[0], y + src_offset[1], z + src_offset[2]);
    }
}

template <typename T,
    size_t X1,
    size_t Y1,
    size_t Z1,
    size_t X2,
    size_t Y2,
    size_t Z2>
INLINE void split_cubo(elem<T> (&dst_p)[X1][Y1][Z1],
    elem<T> (&dst_h)[X1][Y1][Z1],
    const elem2<T> (&src)[X2][Y2][Z2],
    const array<int, 3>& dst_offset,
    const array<int, 3>& src_offset,
    const array<int, 3>& region)
{
    asm volatile("// start split_cubo");

    int gid = get_pid() * THREADS_PER_PE + get_tid();
    int gws = get_maxpid() * THREADS_PER_PE;

    for (int i = gid; i < region[0] * region[1] * region[2]; i += gws) {
        int x = i / region[2] / region[1];
        int y = i / region[2] % region[1];
        int z = i % region[2];

        sync_L1();
        const auto& r = get_elem(src, x + src_offset[0], y + src_offset[1], z + src_offset[2]);

        sync_L1();
        get_elem(dst_p, x + dst_offset[0], y + dst_offset[1], z + dst_offset[2]) = r.p;
        sync_L1();
        get_elem(dst_h, x + dst_offset[0], y + dst_offset[1], z + dst_offset[2]) = r.h;
    }
    asm volatile("// end split_cubo");
    sync_L2();
}

template <typename T,
    size_t X1,
    size_t Y1,
    size_t Z1,
    size_t X2,
    size_t Y2,
    size_t Z2>
INLINE void merge_cubo(elem2<T> (&dst)[X1][Y1][Z1],
    const elem<T> (&src_p)[X2][Y2][Z2],
    const elem<T> (&src_h)[X2][Y2][Z2],
    const array<int, 3>& dst_offset,
    const array<int, 3>& src_offset,
    const array<int, 3>& region)
{
    int gid = get_pid() * THREADS_PER_PE + get_tid();
    int gws = get_maxpid() * THREADS_PER_PE;

    for (int i = gid; i < region[0] * region[1] * region[2]; i += gws) {
        int z = i / region[2] / region[1];
        int y = i / region[2] % region[1];
        int x = i % region[2];

        auto& r = get_elem(dst, x + dst_offset[0], y + dst_offset[1], z + dst_offset[2]);

        r.p = get_elem(src_p, x + src_offset[0], y + src_offset[1], z + src_offset[2]);
        r.h = get_elem(src_h, x + src_offset[0], y + src_offset[1], z + src_offset[2]);
    }
}

template <typename T, size_t X1, size_t Y1, size_t Z1>
INLINE void fill_cubo(T (&dst)[X1][Y1][Z1],
    const T& val,
    const array<int, 3>& dst_offset,
    const array<int, 3>& region)
{
    int gid = get_pid() * THREADS_PER_PE + get_tid();
    int gws = get_maxpid() * THREADS_PER_PE;

    for (int i = gid; i < region[0] * region[1] * region[2]; i += gws) {
        int z = i / region[2] / region[1];
        int y = i / region[2] % region[1];
        int x = i % region[2];

        get_elem(dst, x + dst_offset[0], y + dst_offset[1], z + dst_offset[2]) = val;
    }
}

template <typename T>
using buf = T[NZ + 2 * Ns][NY + 2 * Ns][NX + 2 * Ns];

// constants consts;

template <typename T>
INLINE T d1(T q_m2, T q_m1, T q_0, T q_p1, T q_p2, const constants& consts)
{
    return (q_m2 - consts.num8 * q_m1 + consts.num8 * q_p1 - q_p2) * consts.rev12h1;
}
template <typename T>
INLINE T d2(T q_m2, T q_m1, T q_0, T q_p1, T q_p2, const constants& consts)
{
    return (-q_m2 + consts.num16 * q_m1 - consts.num30 * q_0 + consts.num16 * q_p1 - q_p2) * consts.rev12h2;
}
template <typename T>
INLINE T d3(T q_m2, T q_m1, T q_0, T q_p1, T q_p2, const constants& consts)
{
    return (-q_m2 + consts.num2 * q_m1 - consts.num2 * q_p1 + q_p2) * consts.rev2h3;
}
template <typename T>
INLINE T d4(T q_m2, T q_m1, T q_0, T q_p1, T q_p2, const constants& consts)
{
    return (q_m2 - consts.num4 * q_m1 + consts.num6 * q_0 - consts.num4 * q_p1 + q_p2) * consts.rev1h4;
}

template <typename T>
INLINE T d_x(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d1(q[k][j][i - 2], q[k][j][i - 1], q[k][j][i], q[k][j][i + 1], q[k][j][i + 2], consts);
}
template <typename T>
INLINE T d_y(int i, int j, int k, buf<T>& q, const constants& consts)
{
    return d1(q[k][j - 2][i], q[k][j - 1][i], q[k][j][i], q[k][j + 1][i], q[k][j + 2][i], consts);
}
template <typename T>
INLINE T d_z(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d1(q[k - 2][j][i], q[k - 1][j][i], q[k][j][i], q[k + 1][j][i], q[k + 2][j][i], consts);
}
template <typename T>
INLINE T d_xx(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d2(q[k][j][i - 2], q[k][j][i - 1], q[k][j][i], q[k][j][i + 1], q[k][j][i + 2], consts);
}
template <typename T>
INLINE T d_yy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d2(q[k][j - 2][i], q[k][j - 1][i], q[k][j][i], q[k][j + 1][i], q[k][j + 2][i], consts);
}
template <typename T>
INLINE T d_zz(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d2(q[k - 2][j][i], q[k - 1][j][i], q[k][j][i], q[k + 1][j][i], q[k + 2][j][i], consts);
}
template <typename T>
INLINE T d_xxx(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d3(q[k][j][i - 2], q[k][j][i - 1], q[k][j][i], q[k][j][i + 1], q[k][j][i + 2], consts);
}
template <typename T>
INLINE T d_yyy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d3(q[k][j - 2][i], q[k][j - 1][i], q[k][j][i], q[k][j + 1][i], q[k][j + 2][i], consts);
}
template <typename T>
INLINE T d_zzz(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d3(q[k - 2][j][i], q[k - 1][j][i], q[k][j][i], q[k + 1][j][i], q[k + 2][j][i], consts);
}
template <typename T>
INLINE T d_xxxx(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d4(q[k][j][i - 2], q[k][j][i - 1], q[k][j][i], q[k][j][i + 1], q[k][j][i + 2], consts);
}
template <typename T>
INLINE T d_yyyy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d4(q[k][j - 2][i], q[k][j - 1][i], q[k][j][i], q[k][j + 1][i], q[k][j + 2][i], consts);
}
template <typename T>
INLINE T d_zzzz(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return d4(q[k - 2][j][i], q[k - 1][j][i], q[k][j][i], q[k + 1][j][i], q[k + 2][j][i], consts);
}
template <typename T>
INLINE T d_xy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (q[k + 0][j - 2][i - 2] - consts.num8 * q[k + 0][j - 1][i - 2] + consts.num8 * q[k + 0][j + 1][i - 2] - q[k + 0][j + 2][i - 2] - consts.num8 * q[k + 0][j - 2][i - 1] + consts.num64 * q[k + 0][j - 1][i - 1] - consts.num64 * q[k + 0][j + 1][i - 1] + consts.num8 * q[k + 0][j + 2][i - 1] + consts.num8 * q[k + 0][j - 2][i + 1] - consts.num64 * q[k + 0][j - 1][i + 1] + consts.num64 * q[k + 0][j + 1][i + 1] - consts.num8 * q[k + 0][j + 2][i + 1] - q[k + 0][j - 2][i + 2] + consts.num8 * q[k + 0][j - 1][i + 2] - consts.num8 * q[k + 0][j + 1][i + 2] + q[k + 0][j + 2][i + 2]) * consts.rev144h2;
}
template <typename T>
INLINE T d_xxy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (-q[k + 0][j - 2][i - 2] + consts.num8 * q[k + 0][j - 1][i - 2] - consts.num8 * q[k + 0][j + 1][i - 2] + q[k + 0][j + 2][i - 2] + consts.num16 * q[k + 0][j - 2][i - 1] - consts.num128 * q[k + 0][j - 1][i - 1] + consts.num128 * q[k + 0][j + 1][i - 1] - consts.num16 * q[k + 0][j + 2][i - 1] - consts.num30 * q[k + 0][j - 2][i + 0] + consts.num240 * q[k + 0][j - 1][i + 0] - consts.num240 * q[k + 0][j + 1][i + 0] + consts.num30 * q[k + 0][j + 2][i + 0] + consts.num16 * q[k + 0][j - 2][i + 1] - consts.num128 * q[k + 0][j - 1][i + 1] + consts.num128 * q[k + 0][j + 1][i + 1] - consts.num16 * q[k + 0][j + 2][i + 1] - q[k + 0][j - 2][i + 2] + consts.num8 * q[k + 0][j - 1][i + 2] - consts.num8 * q[k + 0][j + 1][i + 2] + q[k + 0][j + 2][i + 2]) * consts.rev144h3;
}
template <typename T>
INLINE T d_xyy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (-q[k + 0][j - 2][i - 2] + consts.num16 * q[k + 0][j - 1][i - 2] - consts.num30 * q[k + 0][j + 0][i - 2] + consts.num16 * q[k + 0][j + 1][i - 2] - q[k + 0][j + 2][i - 2] + consts.num8 * q[k + 0][j - 2][i - 1] - consts.num128 * q[k + 0][j - 1][i - 1] + consts.num240 * q[k + 0][j + 0][i - 1] - consts.num128 * q[k + 0][j + 1][i - 1] + consts.num8 * q[k + 0][j + 2][i - 1] - consts.num8 * q[k + 0][j - 2][i + 1] + consts.num128 * q[k + 0][j - 1][i + 1] - consts.num240 * q[k + 0][j + 0][i + 1] + consts.num128 * q[k + 0][j + 1][i + 1] - consts.num8 * q[k + 0][j + 2][i + 1] + q[k + 0][j - 2][i + 2] - consts.num16 * q[k + 0][j - 1][i + 2] + consts.num30 * q[k + 0][j + 0][i + 2] - consts.num16 * q[k + 0][j + 1][i + 2] + q[k + 0][j + 2][i + 2]) * consts.rev144h3;
}
template <typename T>
const INLINE T d_xxxy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (-q[k + 0][j - 2][i - 2] + consts.num8 * q[k + 0][j - 1][i - 2] - consts.num8 * q[k + 0][j + 1][i - 2] + q[k + 0][j + 2][i - 2] + consts.num2 * q[k + 0][j - 2][i - 1] - consts.num16 * q[k + 0][j - 1][i - 1] + consts.num16 * q[k + 0][j + 1][i - 1] - consts.num2 * q[k + 0][j + 2][i - 1] - consts.num2 * q[k + 0][j - 2][i + 1] + consts.num16 * q[k + 0][j - 1][i + 1] - consts.num16 * q[k + 0][j + 1][i + 1] + consts.num2 * q[k + 0][j + 2][i + 1] + q[k + 0][j - 2][i + 2] - consts.num8 * q[k + 0][j - 1][i + 2] + consts.num8 * q[k + 0][j + 1][i + 2] - q[k + 0][j + 2][i + 2]) * consts.rev24h4;
}
template <typename T>
INLINE T d_xxyy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (q[k + 0][j - 2][i - 2] - consts.num16 * q[k + 0][j - 1][i - 2] + consts.num30 * q[k + 0][j + 0][i - 2] - consts.num16 * q[k + 0][j + 1][i - 2] + q[k + 0][j + 2][i - 2] - consts.num16 * q[k + 0][j - 2][i - 1] + consts.num256 * q[k + 0][j - 1][i - 1] - consts.num480 * q[k + 0][j + 0][i - 1] + consts.num256 * q[k + 0][j + 1][i - 1] - consts.num16 * q[k + 0][j + 2][i - 1] + consts.num30 * q[k + 0][j - 2][i + 0] - consts.num480 * q[k + 0][j - 1][i + 0] + consts.num900 * q[k + 0][j + 0][i + 0] - consts.num480 * q[k + 0][j + 1][i + 0] + consts.num30 * q[k + 0][j + 2][i + 0] - consts.num16 * q[k + 0][j - 2][i + 1] + consts.num256 * q[k + 0][j - 1][i + 1] - consts.num480 * q[k + 0][j + 0][i + 1] + consts.num256 * q[k + 0][j + 1][i + 1] - consts.num16 * q[k + 0][j + 2][i + 1] + q[k + 0][j - 2][i + 2] - consts.num16 * q[k + 0][j - 1][i + 2] + consts.num30 * q[k + 0][j + 0][i + 2] - consts.num16 * q[k + 0][j + 1][i + 2] + q[k + 0][j + 2][i + 2]) * consts.rev144h4;
}
template <typename T>
INLINE T d_xyyy(int i, int j, int k, const buf<T>& q, const constants& consts)
{
    return (-q[k + 0][j - 2][i - 2] + consts.num2 * q[k + 0][j - 1][i - 2] - consts.num2 * q[k + 0][j + 1][i - 2] + q[k + 0][j + 2][i - 2] + consts.num8 * q[k + 0][j - 2][i - 1] - consts.num16 * q[k + 0][j - 1][i - 1] + consts.num16 * q[k + 0][j + 1][i - 1] - consts.num8 * q[k + 0][j + 2][i - 1] - consts.num8 * q[k + 0][j - 2][i + 1] + consts.num16 * q[k + 0][j - 1][i + 1] - consts.num16 * q[k + 0][j + 1][i + 1] + consts.num8 * q[k + 0][j + 2][i + 1] + q[k + 0][j - 2][i + 2] - consts.num2 * q[k + 0][j - 1][i + 2] + consts.num2 * q[k + 0][j + 1][i + 2] - 1 * q[k + 0][j + 2][i + 2]) * consts.rev24h4;
}

// 差分テーブル
typedef struct
{
    elem<double> x;
    elem<double> y;
    elem<double> xx;
    elem<double> yy;
    elem<double> xy;

    double u_xxx;
    double v_xxx;
    double w_xxx;
    double p_xxx;

    double u_yyy;
    double v_yyy;
    double w_yyy;
    double p_yyy;

    double u_xxy;
    double v_xxy;
    double w_xxy;
    double p_xxy;

    double u_xyy;
    double v_xyy;
    double w_xyy;
    double p_xyy;

    double u_xxxx;
    double v_xxxx;
    double w_xxxx;

    double u_xxxy;
    double v_xxxy;

    double u_xxyy;
    double v_xxyy;
    double w_xxyy;

    double u_xyyy;
    double v_xyyy;

    double u_yyyy;
    double v_yyyy;
    double w_yyyy;
} diffs;

#define assign2(var, suffix, zz, mem1, mem2)            \
    {                                                   \
        auto e = d_##suffix(ix, iy, zz, buf.p, consts); \
        var.mem1##_##suffix = e.mem1;                   \
        var.mem2##_##suffix = e.mem2;                   \
    }
#define assign3(var, suffix, zz, mem1, mem2, mem3)      \
    {                                                   \
        auto e = d_##suffix(ix, iy, zz, buf.p, consts); \
        var.mem1##_##suffix = e.mem1;                   \
        var.mem2##_##suffix = e.mem2;                   \
        var.mem3##_##suffix = e.mem3;                   \
    }
#define assign4(var, suffix, zz, mem1, mem2, mem3, mem4) \
    {                                                    \
        auto e = d_##suffix(ix, iy, zz, buf.p, consts);  \
        var.mem1##_##suffix = e.mem1;                    \
        var.mem2##_##suffix = e.mem2;                    \
        var.mem3##_##suffix = e.mem3;                    \
        var.mem4##_##suffix = e.mem4;                    \
    }

#define def2(expr, suffix, mem1, mem2) \
    double mem1##_##suffix;            \
    double mem2##_##suffix;            \
    {                                  \
        auto e = expr;                 \
        mem1##_##suffix = e.mem1;      \
        mem2##_##suffix = e.mem2;      \
    }
#define def3(expr, suffix, mem1, mem2, mem3) \
    double mem1##_##suffix;                  \
    double mem2##_##suffix;                  \
    double mem3##_##suffix;                  \
    {                                        \
        auto e = expr;                       \
        mem1##_##suffix = e.mem1;            \
        mem2##_##suffix = e.mem2;            \
        mem3##_##suffix = e.mem3;            \
    }
#define def4(expr, suffix, mem1, mem2, mem3, mem4) \
    double mem1##_##suffix;                        \
    double mem2##_##suffix;                        \
    double mem3##_##suffix;                        \
    double mem4##_##suffix;                        \
    {                                              \
        auto e = expr;                             \
        mem1##_##suffix = e.mem1;                  \
        mem2##_##suffix = e.mem2;                  \
        mem3##_##suffix = e.mem3;                  \
        mem4##_##suffix = e.mem4;                  \
    }

// 1ステップ更新
INLINE void step(
    int ix, int iy, int bz, int WZ,
    local_buf<double>& buf,
    local_buf<double>& res)
{
    TakeProfile<false, true> prof_all_(auxPerf[1]);

    static_assert(sizeof(diffs[2 * Ns + 1]) * THREADS_PER_PE <= 10 * 1024,
        "Local memory usage is too large");

    struct local_memory {
        diffs ds[THREADS_PER_PE][2 * Ns + 1];
        constants consts;
    };

    local_memory& lm = *(local_memory*)(20 * 1024 - sizeof(local_memory));
    auto& ds = lm.ds[get_tid()];
    auto& consts = lm.consts;

#pragma nounroll
    for (int k = 0; k < 2 * Ns; k++) {
        sync_L2();

        asm volatile("// precalc start");

        ds[k].x = d_x(ix, iy, bz * WZ + k, buf.p, consts);
        ds[k].y = d_y(ix, iy, bz * WZ + k, buf.p, consts);
        ds[k].xx = d_xx(ix, iy, bz * WZ + k, buf.p, consts);
        ds[k].xy = d_xy(ix, iy, bz * WZ + k, buf.p, consts);
        ds[k].yy = d_yy(ix, iy, bz * WZ + k, buf.p, consts);

        assign4(ds[k], xxx, bz * WZ + k, u, v, w, p);
        assign4(ds[k], xxy, bz * WZ + k, u, v, w, p);
        assign4(ds[k], xyy, bz * WZ + k, u, v, w, p);
        assign4(ds[k], yyy, bz * WZ + k, u, v, w, p);
        assign3(ds[k], xxxx, bz * WZ + k, u, v, w);
        assign2(ds[k], xxxy, bz * WZ + k, u, v);
        assign3(ds[k], xxyy, bz * WZ + k, u, v, w);
        assign2(ds[k], xyyy, bz * WZ + k, u, v);
        assign3(ds[k], yyyy, bz * WZ + k, u, v, w);

        asm volatile("// precalc end");
    }

#pragma nounroll
    for (int iz = Ns + bz * WZ; iz < Ns + bz * WZ + WZ; iz++) {
        sync_L2();

        asm volatile("// step start");

#define shift(x) ((x) % (2 * Ns + 1))

        int k_p2 = shift(iz + 2);
        int k_p1 = shift(iz + 1);
        int k_0 = shift(iz);
        int k_m1 = shift(iz - 1);
        int k_m2 = shift(iz - 2);

        ds[k_p2].x = d_x(ix, iy, iz + 2, buf.p, consts);
        ds[k_p2].y = d_y(ix, iy, iz + 2, buf.p, consts);
        ds[k_p2].xx = d_xx(ix, iy, iz + 2, buf.p, consts);
        ds[k_p2].yy = d_yy(ix, iy, iz + 2, buf.p, consts);
        ds[k_p2].xy = d_xy(ix, iy, iz + 2, buf.p, consts);
        assign4(ds[k_p2], xxx, iz + 2, u, v, w, p);
        assign4(ds[k_p2], xxy, iz + 2, u, v, w, p);
        assign4(ds[k_p2], xyy, iz + 2, u, v, w, p);
        assign4(ds[k_p2], yyy, iz + 2, u, v, w, p);
        assign3(ds[k_p2], xxxx, iz + 2, u, v, w);
        assign2(ds[k_p2], xxxy, iz + 2, u, v);
        assign3(ds[k_p2], xxyy, iz + 2, u, v, w);
        assign2(ds[k_p2], xyyy, iz + 2, u, v);
        assign3(ds[k_p2], yyyy, iz + 2, u, v, w);

        auto& rp = buf.p[iz][iy][ix];
        double r = rp.r;
        double u = rp.u;
        double v = rp.v;
        double w = rp.w;
        double p = rp.p;

        double r_x = ds[k_0].x.r;
        double u_x = ds[k_0].x.u;
        double v_x = ds[k_0].x.v;
        double w_x = ds[k_0].x.w;
        double p_x = ds[k_0].x.p;
        double r_y = ds[k_0].y.r;
        double u_y = ds[k_0].y.u;
        double v_y = ds[k_0].y.v;
        double w_y = ds[k_0].y.w;
        double p_y = ds[k_0].y.p;

        auto e1 = d_z(ix, iy, iz, buf.p, consts);
        double r_z = e1.r;
        double u_z = e1.u;
        double v_z = e1.v;
        double w_z = e1.w;
        double p_z = e1.p;

        double r_xx = ds[k_0].xx.r;
        double u_xx = ds[k_0].xx.u;
        double v_xx = ds[k_0].xx.v;
        double w_xx = ds[k_0].xx.w;
        double p_xx = ds[k_0].xx.p;
        double r_yy = ds[k_0].yy.r;
        double u_yy = ds[k_0].yy.u;
        double v_yy = ds[k_0].yy.v;
        double w_yy = ds[k_0].yy.w;
        double p_yy = ds[k_0].yy.p;
        auto e2 = d_zz(ix, iy, iz, buf.p, consts);
        double r_zz = e2.r;
        double u_zz = e2.u;
        double v_zz = e2.v;
        double w_zz = e2.w;
        double p_zz = e2.p;
        double r_xy = ds[k_0].xy.r;
        double u_xy = ds[k_0].xy.u;
        double v_xy = ds[k_0].xy.v;
        double w_xy = ds[k_0].xy.w;
        double p_xy = ds[k_0].xy.p;

        auto e3 = d1(ds[k_m2].y, ds[k_m1].y, ds[k_0].y, ds[k_p1].y, ds[k_p2].y, consts);
        double r_yz = e3.r;
        double u_yz = e3.u;
        double v_yz = e3.v;
        double w_yz = e3.w;
        double p_yz = e3.p;

        auto e4 = d1(ds[k_m2].x, ds[k_m1].x, ds[k_0].x, ds[k_p1].x, ds[k_p2].x, consts);
        double r_xz = e4.r;
        double u_xz = e4.u;
        double v_xz = e4.v;
        double w_xz = e4.w;
        double p_xz = e4.p;

        double u_xxx = ds[k_0].u_xxx;
        double v_xxx = ds[k_0].v_xxx;
        double w_xxx = ds[k_0].w_xxx;
        double p_xxx = ds[k_0].p_xxx;
        double u_yyy = ds[k_0].u_yyy;
        double v_yyy = ds[k_0].v_yyy;
        double w_yyy = ds[k_0].w_yyy;
        double p_yyy = ds[k_0].p_yyy;
        def4((d_zzz(ix, iy, iz, buf.p, consts)), zzz, u, v, w, p);
        double u_xxy = ds[k_0].u_xxy;
        double v_xxy = ds[k_0].v_xxy;
        double w_xxy = ds[k_0].w_xxy;
        double p_xxy = ds[k_0].p_xxy;
        def4((d1(ds[k_m2].xx, ds[k_m1].xx, ds[k_0].xx, ds[k_p1].xx, ds[k_p2].xx, consts)), xxz, u, v, w, p);
        double u_xyy = ds[k_0].u_xyy;
        double v_xyy = ds[k_0].v_xyy;
        double w_xyy = ds[k_0].w_xyy;
        double p_xyy = ds[k_0].p_xyy;
        def4((d2(ds[k_m2].x, ds[k_m1].x, ds[k_0].x, ds[k_p1].x, ds[k_p2].x, consts)), xzz, u, v, w, p);
        def4((d1(ds[k_m2].yy, ds[k_m1].yy, ds[k_0].yy, ds[k_p1].yy, ds[k_p2].yy, consts)), yyz, u, v, w, p);
        def4((d2(ds[k_m2].y, ds[k_m1].y, ds[k_0].y, ds[k_p1].y, ds[k_p2].y, consts)), yzz, u, v, w, p);
        def3((d1(ds[k_m2].xy, ds[k_m1].xy, ds[k_0].xy, ds[k_p1].xy, ds[k_p2].xy, consts)), xyz, u, v, w);

        double u_xxxx = ds[k_0].u_xxxx;
        double v_xxxx = ds[k_0].v_xxxx;
        double w_xxxx = ds[k_0].w_xxxx;

        double u_xxxy = ds[k_0].u_xxxy;
        double v_xxxy = ds[k_0].v_xxxy;

        double u_xxyy = ds[k_0].u_xxyy;
        double v_xxyy = ds[k_0].v_xxyy;
        double w_xxyy = ds[k_0].w_xxyy;

        double u_xyyy = ds[k_0].u_xyyy;
        double v_xyyy = ds[k_0].v_xyyy;

        double u_yyyy = ds[k_0].u_yyyy;
        double v_yyyy = ds[k_0].v_yyyy;
        double w_yyyy = ds[k_0].w_yyyy;

        double u_xxxz = d1(ds[k_m2].u_xxx, ds[k_m1].u_xxx, ds[k_0].u_xxx, ds[k_p1].u_xxx, ds[k_p2].u_xxx, consts);
        double w_xxxz = d1(ds[k_m2].w_xxx, ds[k_m1].w_xxx, ds[k_0].w_xxx, ds[k_p1].w_xxx, ds[k_p2].w_xxx, consts);

        double v_xxyz = d1(ds[k_m2].v_xxy, ds[k_m1].v_xxy, ds[k_0].v_xxy, ds[k_p1].v_xxy, ds[k_p2].v_xxy, consts);
        double w_xxyz = d1(ds[k_m2].w_xxy, ds[k_m1].w_xxy, ds[k_0].w_xxy, ds[k_p1].w_xxy, ds[k_p2].w_xxy, consts);

        double u_xyyz = d1(ds[k_m2].u_xyy, ds[k_m1].u_xyy, ds[k_0].u_xyy, ds[k_p1].u_xyy, ds[k_p2].u_xyy, consts);
        double w_xyyz = d1(ds[k_m2].w_xyy, ds[k_m1].w_xyy, ds[k_0].w_xyy, ds[k_p1].w_xyy, ds[k_p2].w_xyy, consts);

        def3((d2(ds[k_m2].xx, ds[k_m1].xx, ds[k_0].xx, ds[k_p1].xx, ds[k_p2].xx, consts)), xxzz, u, v, w);
        def2((d2(ds[k_m2].xy, ds[k_m1].xy, ds[k_0].xy, ds[k_p1].xy, ds[k_p2].xy, consts)), xyzz, u, v);
        def2((d2(ds[k_m2].x, ds[k_m1].x, ds[k_0].x, ds[k_p1].x, ds[k_p2].x, consts)), xzzz, u, w);

        double v_yyyz = d1(ds[k_m2].v_yyy, ds[k_m1].v_yyy, ds[k_0].v_yyy, ds[k_p1].v_yyy, ds[k_p2].v_yyy, consts);
        double w_yyyz = d1(ds[k_m2].w_yyy, ds[k_m1].w_yyy, ds[k_0].w_yyy, ds[k_p1].w_yyy, ds[k_p2].w_yyy, consts);

        def3((d2(ds[k_m2].yy, ds[k_m1].yy, ds[k_0].yy, ds[k_p1].yy, ds[k_p2].yy, consts)), yyzz, u, v, w);
        def2((d2(ds[k_m2].y, ds[k_m1].y, ds[k_0].y, ds[k_p1].y, ds[k_p2].y, consts)), yzzz, v, w);
        def3((d_zzzz(ix, iy, iz, buf.p, consts)), zzzz, u, v, w);

        if (get_pid() >= 9999) {
            sync();
        }

        //-----

        double vis1 = consts.num4 / consts.num3 * u_xx + consts.num1 / consts.num3 * v_xy + consts.num1 / consts.num3 * w_xz + u_yy + u_zz;
        double vis2 = consts.num1 / consts.num3 * u_xy + consts.num4 / consts.num3 * v_yy + consts.num1 / consts.num3 * w_yz + v_xx + v_zz;
        double vis3 = consts.num1 / consts.num3 * u_xz + consts.num1 / consts.num3 * v_yz + consts.num4 / consts.num3 * w_zz + w_xx + w_yy;
        double vis1_x = consts.num4 / consts.num3 * u_xxx + consts.num1 / consts.num3 * v_xxy + consts.num1 / consts.num3 * w_xxz + u_xyy + u_xzz;
        double vis2_x = consts.num1 / consts.num3 * u_xxy + consts.num4 / consts.num3 * v_xyy + consts.num1 / consts.num3 * w_xyz + v_xzz + v_xxx;
        double vis3_x = consts.num1 / consts.num3 * u_xxz + consts.num1 / consts.num3 * v_xyz + consts.num4 / consts.num3 * w_xzz + w_xyy + w_xxx;
        double vis1_y = consts.num4 / consts.num3 * u_xxy + consts.num1 / consts.num3 * v_xyy + consts.num1 / consts.num3 * w_xyz + u_yzz + u_yyy;
        double vis2_y = consts.num1 / consts.num3 * u_xyy + consts.num4 / consts.num3 * v_yyy + consts.num1 / consts.num3 * w_yyz + v_xxy + v_yzz;
        double vis3_y = consts.num1 / consts.num3 * u_xyz + consts.num1 / consts.num3 * v_yyz + consts.num4 / consts.num3 * w_yzz + w_xxy + w_yyy;
        double vis1_z = consts.num4 / consts.num3 * u_xxz + consts.num1 / consts.num3 * v_xyz + consts.num1 / consts.num3 * w_xzz + u_yyz + u_zzz;
        double vis2_z = consts.num1 / consts.num3 * u_xyz + consts.num4 / consts.num3 * v_yyz + consts.num1 / consts.num3 * w_yzz + v_xxz + v_zzz;
        double vis3_z = consts.num1 / consts.num3 * u_xzz + consts.num1 / consts.num3 * v_yzz + consts.num4 / consts.num3 * w_zzz + w_xxz + w_yyz;
        double vis1_xx = consts.num4 / consts.num3 * u_xxxx + consts.num1 / consts.num3 * v_xxxy + consts.num1 / consts.num3 * w_xxxz + u_xxyy + u_xxzz;
        double vis2_xx = consts.num1 / consts.num3 * u_xxxy + consts.num4 / consts.num3 * v_xxyy + consts.num1 / consts.num3 * w_xxyz + v_xxzz + v_xxxx;
        double vis3_xx = consts.num1 / consts.num3 * u_xxxz + consts.num1 / consts.num3 * v_xxyz + consts.num4 / consts.num3 * w_xxzz + w_xxyy + w_xxxx;
        double vis1_yy = consts.num4 / consts.num3 * u_xxyy + consts.num1 / consts.num3 * v_xyyy + consts.num1 / consts.num3 * w_xyyz + u_yyzz + u_yyyy;
        double vis2_yy = consts.num1 / consts.num3 * u_xyyy + consts.num4 / consts.num3 * v_yyyy + consts.num1 / consts.num3 * w_yyyz + v_xxyy + v_yyzz;
        double vis3_yy = consts.num1 / consts.num3 * u_xyyz + consts.num1 / consts.num3 * v_yyyz + consts.num4 / consts.num3 * w_yyzz + w_xxyy + w_yyyy;
        double vis1_zz = consts.num4 / consts.num3 * u_xxzz + consts.num1 / consts.num3 * v_xyzz + consts.num1 / consts.num3 * w_xzzz + u_yyzz + u_zzzz;
        double vis2_zz = consts.num1 / consts.num3 * u_xyzz + consts.num4 / consts.num3 * v_yyzz + consts.num1 / consts.num3 * w_yzzz + v_xxzz + v_zzzz;
        double vis3_zz = consts.num1 / consts.num3 * u_xzzz + consts.num1 / consts.num3 * v_yzzz + consts.num4 / consts.num3 * w_zzzz + w_xxzz + w_yyzz;
        double vis1_xy = consts.num4 / consts.num3 * u_xxxy + consts.num1 / consts.num3 * v_xxyy + consts.num1 / consts.num3 * w_xxyz + u_xyzz + u_xyyy;
        double vis2_xy = consts.num1 / consts.num3 * u_xxyy + consts.num4 / consts.num3 * v_xyyy + consts.num1 / consts.num3 * w_xyyz + v_xyzz + v_xxxy;
        double vis2_yz = consts.num1 / consts.num3 * u_xyyz + consts.num4 / consts.num3 * v_yyyz + consts.num1 / consts.num3 * w_yyzz + v_xxyz + v_yzzz;
        double vis3_yz = consts.num1 / consts.num3 * u_xyzz + consts.num1 / consts.num3 * v_yyzz + consts.num4 / consts.num3 * w_yzzz + w_xxyz + w_yyyz;
        double vis1_xz = consts.num4 / consts.num3 * u_xxxz + consts.num1 / consts.num3 * v_xxyz + consts.num1 / consts.num3 * w_xxzz + u_xyyz + u_xzzz;
        double vis3_xz = consts.num1 / consts.num3 * u_xxzz + consts.num1 / consts.num3 * v_xyzz + consts.num4 / consts.num3 * w_xzzz + w_xyyz + w_xxxz;

        double rr1 = consts.num1 / r;
        double rr2 = consts.num1 / (r * r);
        double rr3 = consts.num1 / (r * r * r);

        double r_t = -r * u_x - r * v_y - r * w_z - r_x * u - r_y * v - r_z * w;
        double u_t = consts.c * rr1 * vis1 - p_x * rr1 - u * u_x - u_y * v - u_z * w;
        double v_t = consts.c * rr1 * vis2 - p_y * rr1 - u * v_x - v * v_y - v_z * w;
        double w_t = consts.c * rr1 * vis3 - p_z * rr1 - u * w_x - v * w_y - w * w_z;
        double p_t = -consts.c2 * u * vis1 - consts.c2 * v * vis2 - consts.c2 * vis3 * w - consts.gm * p * u_x - consts.gm * p * v_y - consts.gm * p * w_z - p_x * u - p_y * v - p_z * w;

        double r_tx = -r * u_xx - r * v_xy - r * w_xz - consts.num2 * r_x * u_x - r_x * v_y - r_x * w_z - r_xy * v - r_xz * w - r_xx * u - r_y * v_x - r_z * w_x;
        double u_tx = -consts.c * rr2 * r_x * vis1 + consts.c * rr1 * vis1_x + p_x * rr2 * r_x - p_xx * rr1 - u * u_xx - (u_x * u_x) - u_xy * v - u_xz * w - u_y * v_x - u_z * w_x;
        double v_tx = -consts.c * rr2 * r_x * vis2 + consts.c * rr1 * vis2_x - p_xy * rr1 + p_y * rr2 * r_x - u * v_xx - u_x * v_x - v * v_xy - v_x * v_y - v_xz * w - v_z * w_x;
        double w_tx = -consts.c * rr2 * r_x * vis3 + consts.c * rr1 * vis3_x - p_xz * rr1 + p_z * rr2 * r_x - u * w_xx - u_x * w_x - v * w_xy - v_x * w_y - w * w_xz - w_x * w_z;
        double p_tx = -consts.c2 * u * vis1_x - consts.c2 * u_x * vis1 - consts.c2 * v * vis2_x - consts.c2 * v_x * vis2 - consts.c2 * vis3 * w_x - consts.c2 * vis3_x * w - consts.gm * p * u_xx - consts.gm * p * v_xy - consts.gm * p * w_xz - consts.gm * p_x * u_x - consts.gm * p_x * v_y - consts.gm * p_x * w_z - p_x * u_x - p_xy * v - p_xz * w - p_xx * u - p_y * v_x - p_z * w_x;

        double r_ty = -r * u_xy - r * v_yy - r * w_yz - r_x * u_y - r_xy * u - r_y * u_x - consts.num2 * r_y * v_y - r_y * w_z - r_yz * w - r_yy * v - r_z * w_y;
        double u_ty = -consts.c * rr2 * r_y * vis1 + consts.c * rr1 * vis1_y + p_x * rr2 * r_y - p_xy * rr1 - u * u_xy - u_x * u_y - u_y * v_y - u_yz * w - u_yy * v - u_z * w_y;
        double v_ty = -consts.c * rr2 * r_y * vis2 + consts.c * rr1 * vis2_y + p_y * rr2 * r_y - p_yy * rr1 - u * v_xy - u_y * v_x - v * v_yy - (v_y * v_y) - v_yz * w - v_z * w_y;
        double w_ty = -consts.c * rr2 * r_y * vis3 + consts.c * rr1 * vis3_y - p_yz * rr1 + p_z * rr2 * r_y - u * w_xy - u_y * w_x - v * w_yy - v_y * w_y - w * w_yz - w_y * w_z;
        double p_ty = -consts.c2 * u * vis1_y - consts.c2 * u_y * vis1 - consts.c2 * v * vis2_y - consts.c2 * v_y * vis2 - consts.c2 * vis3 * w_y - consts.c2 * vis3_y * w - consts.gm * p * u_xy - consts.gm * p * v_yy - consts.gm * p * w_yz - consts.gm * p_y * u_x - consts.gm * p_y * v_y - consts.gm * p_y * w_z - p_x * u_y - p_xy * u - p_y * v_y - p_yz * w - p_yy * v - p_z * w_y;

        double r_tz = -r * u_xz - r * v_yz - r * w_zz - r_x * u_z - r_xz * u - r_y * v_z - r_yz * v - r_z * u_x - r_z * v_y - consts.num2 * r_z * w_z - r_zz * w;
        double u_tz = -consts.c * rr2 * r_z * vis1 + consts.c * rr1 * vis1_z + p_x * rr2 * r_z - p_xz * rr1 - u * u_xz - u_x * u_z - u_y * v_z - u_yz * v - u_z * w_z - u_zz * w;
        double v_tz = -consts.c * rr2 * r_z * vis2 + consts.c * rr1 * vis2_z + p_y * rr2 * r_z - p_yz * rr1 - u * v_xz - u_z * v_x - v * v_yz - v_y * v_z - v_z * w_z - v_zz * w;
        double w_tz = -consts.c * rr2 * r_z * vis3 + consts.c * rr1 * vis3_z + p_z * rr2 * r_z - p_zz * rr1 - u * w_xz - u_z * w_x - v * w_yz - v_z * w_y - w * w_zz - (w_z * w_z);
        double p_tz = -consts.c2 * u * vis1_z - consts.c2 * u_z * vis1 - consts.c2 * v * vis2_z - consts.c2 * v_z * vis2 - consts.c2 * vis3 * w_z - consts.c2 * vis3_z * w - consts.gm * p * u_xz - consts.gm * p * v_yz - consts.gm * p * w_zz - consts.gm * p_z * u_x - consts.gm * p_z * v_y - consts.gm * p_z * w_z - p_x * u_z - p_xz * u - p_y * v_z - p_yz * v - p_z * w_z - p_zz * w;

        double u_txx = consts.num2 * consts.c * rr3 * (r_x * r_x) * vis1 - consts.num2 * consts.c * rr2 * r_x * vis1_x - consts.c * rr2 * r_xx * vis1 + consts.c * rr1 * vis1_xx - consts.num2 * p_x * rr3 * (r_x * r_x) + p_x * rr2 * r_xx + consts.num2 * p_xx * rr2 * r_x - p_xxx * rr1 - u * u_xxx - consts.num3 * u_x * u_xx - consts.num2 * u_xy * v_x - consts.num2 * u_xz * w_x - u_xxy * v - u_xxz * w - u_y * v_xx - u_z * w_xx;
        double v_txx = consts.num2 * consts.c * rr3 * (r_x * r_x) * vis2 - consts.num2 * consts.c * rr2 * r_x * vis2_x - consts.c * rr2 * r_xx * vis2 + consts.c * rr1 * vis2_xx + consts.num2 * p_xy * rr2 * r_x - p_xxy * rr1 - consts.num2 * p_y * rr3 * (r_x * r_x) + p_y * rr2 * r_xx - u * v_xxx - consts.num2 * u_x * v_xx - u_xx * v_x - v * v_xxy - consts.num2 * v_x * v_xy - consts.num2 * v_xz * w_x - v_xx * v_y - v_xxz * w - v_z * w_xx;
        double w_txx = consts.num2 * consts.c * rr3 * (r_x * r_x) * vis3 - consts.num2 * consts.c * rr2 * r_x * vis3_x - consts.c * rr2 * r_xx * vis3 + consts.c * rr1 * vis3_xx + consts.num2 * p_xz * rr2 * r_x - p_xxz * rr1 - consts.num2 * p_z * rr3 * (r_x * r_x) + p_z * rr2 * r_xx - u * w_xxx - consts.num2 * u_x * w_xx - u_xx * w_x - v * w_xxy - consts.num2 * v_x * w_xy - v_xx * w_y - w * w_xxz - consts.num2 * w_x * w_xz - w_xx * w_z;

        double u_tyy = consts.num2 * consts.c * rr3 * (r_y * r_y) * vis1 - consts.num2 * consts.c * rr2 * r_y * vis1_y - consts.c * rr2 * r_yy * vis1 + consts.c * rr1 * vis1_yy - consts.num2 * p_x * rr3 * (r_y * r_y) + p_x * rr2 * r_yy + consts.num2 * p_xy * rr2 * r_y - p_xyy * rr1 - u * u_xyy - u_x * u_yy - consts.num2 * u_xy * u_y - u_y * v_yy - consts.num2 * u_yz * w_y - consts.num2 * u_yy * v_y - u_yyz * w - u_yyy * v - u_z * w_yy;
        double v_tyy = consts.num2 * consts.c * rr3 * (r_y * r_y) * vis2 - consts.num2 * consts.c * rr2 * r_y * vis2_y - consts.c * rr2 * r_yy * vis2 + consts.c * rr1 * vis2_yy - consts.num2 * p_y * rr3 * (r_y * r_y) + p_y * rr2 * r_yy + consts.num2 * p_yy * rr2 * r_y - p_yyy * rr1 - u * v_xyy - consts.num2 * u_y * v_xy - u_yy * v_x - v * v_yyy - consts.num3 * v_y * v_yy - consts.num2 * v_yz * w_y - v_yyz * w - v_z * w_yy;
        double w_tyy = consts.num2 * consts.c * rr3 * (r_y * r_y) * vis3 - consts.num2 * consts.c * rr2 * r_y * vis3_y - consts.c * rr2 * r_yy * vis3 + consts.c * rr1 * vis3_yy + consts.num2 * p_yz * rr2 * r_y - p_yyz * rr1 - consts.num2 * p_z * rr3 * (r_y * r_y) + p_z * rr2 * r_yy - u * w_xyy - consts.num2 * u_y * w_xy - u_yy * w_x - v * w_yyy - consts.num2 * v_y * w_yy - v_yy * w_y - w * w_yyz - consts.num2 * w_y * w_yz - w_yy * w_z;

        double u_tzz = consts.num2 * consts.c * rr3 * (r_z * r_z) * vis1 - consts.num2 * consts.c * rr2 * r_z * vis1_z - consts.c * rr2 * r_zz * vis1 + consts.c * rr1 * vis1_zz - consts.num2 * p_x * rr3 * (r_z * r_z) + p_x * rr2 * r_zz + consts.num2 * p_xz * rr2 * r_z - p_xzz * rr1 - u * u_xzz - u_x * u_zz - consts.num2 * u_xz * u_z - u_y * v_zz - consts.num2 * u_yz * v_z - u_yzz * v - u_z * w_zz - consts.num2 * u_zz * w_z - u_zzz * w;
        double v_tzz = consts.num2 * consts.c * rr3 * (r_z * r_z) * vis2 - consts.num2 * consts.c * rr2 * r_z * vis2_z - consts.c * rr2 * r_zz * vis2 + consts.c * rr1 * vis2_zz - consts.num2 * p_y * rr3 * (r_z * r_z) + p_y * rr2 * r_zz + consts.num2 * p_yz * rr2 * r_z - p_yzz * rr1 - u * v_xzz - consts.num2 * u_z * v_xz - u_zz * v_x - v * v_yzz - v_y * v_zz - consts.num2 * v_yz * v_z - v_z * w_zz - consts.num2 * v_zz * w_z - v_zzz * w;
        double w_tzz = consts.num2 * consts.c * rr3 * (r_z * r_z) * vis3 - consts.num2 * consts.c * rr2 * r_z * vis3_z - consts.c * rr2 * r_zz * vis3 + consts.c * rr1 * vis3_zz - consts.num2 * p_z * rr3 * (r_z * r_z) + p_z * rr2 * r_zz + consts.num2 * p_zz * rr2 * r_z - p_zzz * rr1 - u * w_xzz - consts.num2 * u_z * w_xz - u_zz * w_x - v * w_yzz - consts.num2 * v_z * w_yz - v_zz * w_y - w * w_zzz - consts.num3 * w_z * w_zz;

        double u_txy = consts.num2 * consts.c * rr3 * r_x * r_y * vis1 - consts.c * rr2 * r_x * vis1_y - consts.c * rr2 * r_xy * vis1 - consts.c * rr2 * r_y * vis1_x + consts.c * rr1 * vis1_xy - consts.num2 * p_x * rr3 * r_x * r_y + p_x * rr2 * r_xy + p_xy * rr2 * r_x + p_xx * rr2 * r_y - p_xxy * rr1 - u * u_xxy - consts.num2 * u_x * u_xy - u_xy * v_y - u_xyz * w - u_xyy * v - u_xz * w_y - u_xx * u_y - u_y * v_xy - u_yz * w_x - u_yy * v_x - u_z * w_xy;
        double v_txy = consts.num2 * consts.c * rr3 * r_x * r_y * vis2 - consts.c * rr2 * r_x * vis2_y - consts.c * rr2 * r_xy * vis2 - consts.c * rr2 * r_y * vis2_x + consts.c * rr1 * vis2_xy + p_xy * rr2 * r_y - p_xyy * rr1 - consts.num2 * p_y * rr3 * r_x * r_y + p_y * rr2 * r_xy + p_yy * rr2 * r_x - u * v_xxy - u_x * v_xy - u_xy * v_x - u_y * v_xx - v * v_xyy - v_x * v_yy - consts.num2 * v_xy * v_y - v_xyz * w - v_xz * w_y - v_yz * w_x - v_z * w_xy;

        double v_tyz = consts.num2 * consts.c * rr3 * r_y * r_z * vis2 - consts.c * rr2 * r_y * vis2_z - consts.c * rr2 * r_yz * vis2 - consts.c * rr2 * r_z * vis2_y + consts.c * rr1 * vis2_yz - consts.num2 * p_y * rr3 * r_y * r_z + p_y * rr2 * r_yz + p_yz * rr2 * r_y + p_yy * rr2 * r_z - p_yyz * rr1 - u * v_xyz - u_y * v_xz - u_yz * v_x - u_z * v_xy - v * v_yyz - consts.num2 * v_y * v_yz - v_yz * w_z - v_yzz * w - v_yy * v_z - v_z * w_yz - v_zz * w_y;
        double w_tyz = consts.num2 * consts.c * rr3 * r_y * r_z * vis3 - consts.c * rr2 * r_y * vis3_z - consts.c * rr2 * r_yz * vis3 - consts.c * rr2 * r_z * vis3_y + consts.c * rr1 * vis3_yz + p_yz * rr2 * r_z - p_yzz * rr1 - consts.num2 * p_z * rr3 * r_y * r_z + p_z * rr2 * r_yz + p_zz * rr2 * r_y - u * w_xyz - u_y * w_xz - u_yz * w_x - u_z * w_xy - v * w_yyz - v_y * w_yz - v_yz * w_y - v_z * w_yy - w * w_yzz - w_y * w_zz - consts.num2 * w_yz * w_z;

        double u_txz = consts.num2 * consts.c * rr3 * r_x * r_z * vis1 - consts.c * rr2 * r_x * vis1_z - consts.c * rr2 * r_xz * vis1 - consts.c * rr2 * r_z * vis1_x + consts.c * rr1 * vis1_xz - consts.num2 * p_x * rr3 * r_x * r_z + p_x * rr2 * r_xz + p_xz * rr2 * r_x + p_xx * rr2 * r_z - p_xxz * rr1 - u * u_xxz - consts.num2 * u_x * u_xz - u_xy * v_z - u_xyz * v - u_xz * w_z - u_xzz * w - u_xx * u_z - u_y * v_xz - u_yz * v_x - u_z * w_xz - u_zz * w_x;
        double w_txz = consts.num2 * consts.c * rr3 * r_x * r_z * vis3 - consts.c * rr2 * r_x * vis3_z - consts.c * rr2 * r_xz * vis3 - consts.c * rr2 * r_z * vis3_x + consts.c * rr1 * vis3_xz + p_xz * rr2 * r_z - p_xzz * rr1 - consts.num2 * p_z * rr3 * r_x * r_z + p_z * rr2 * r_xz + p_zz * rr2 * r_x - u * w_xxz - u_x * w_xz - u_xz * w_x - u_z * w_xx - v * w_xyz - v_x * w_yz - v_xz * w_y - v_z * w_xy - w * w_xzz - w_x * w_zz - consts.num2 * w_xz * w_z;

        double vis1_t = consts.num4 / consts.num3 * u_txx + consts.num1 / consts.num3 * v_txy + consts.num1 / consts.num3 * w_txz + u_tyy + u_tzz;
        double vis2_t = consts.num1 / consts.num3 * u_txy + consts.num4 / consts.num3 * v_tyy + consts.num1 / consts.num3 * w_tyz + v_txx + v_tzz;
        double vis3_t = consts.num1 / consts.num3 * u_txz + consts.num1 / consts.num3 * v_tyz + consts.num4 / consts.num3 * w_tzz + w_txx + w_tyy;

        double r_tt = -r * u_tx - r * v_ty - r * w_tz - r_t * u_x - r_t * v_y - r_t * w_z - r_tx * u - r_ty * v - r_tz * w - r_x * u_t - r_y * v_t - r_z * w_t;
        double u_tt = -consts.c * rr2 * r_t * vis1 + consts.c * rr1 * vis1_t - p_tx * rr1 + p_x * rr2 * r_t - u * u_tx - u_t * u_x - u_ty * v - u_tz * w - u_y * v_t - u_z * w_t;
        double v_tt = -consts.c * rr2 * r_t * vis2 + consts.c * rr1 * vis2_t - p_ty * rr1 + p_y * rr2 * r_t - u * v_tx - u_t * v_x - v * v_ty - v_t * v_y - v_tz * w - v_z * w_t;
        double w_tt = -consts.c * rr2 * r_t * vis3 + consts.c * rr1 * vis3_t - p_tz * rr1 + p_z * rr2 * r_t - u * w_tx - u_t * w_x - v * w_ty - v_t * w_y - w * w_tz - w_t * w_z;
        double p_tt = -consts.c2 * u * vis1_t - consts.c2 * u_t * vis1 - consts.c2 * v * vis2_t - consts.c2 * v_t * vis2 - consts.c2 * vis3 * w_t - consts.c2 * vis3_t * w - consts.gm * p * u_tx - consts.gm * p * v_ty - consts.gm * p * w_tz - consts.gm * p_t * u_x - consts.gm * p_t * v_y - consts.gm * p_t * w_z - p_tx * u - p_ty * v - p_tz * w - p_x * u_t - p_y * v_t - p_z * w_t;

        get_elem(res.p, ix, iy, iz - Ns)
            = get_elem(buf.h, ix, iy, iz) + elem<double>{ r_t, u_t, v_t, w_t, p_t } * consts.num3 * consts.dt * (consts.num1 / consts.num2) + consts.num5 * consts.dt * consts.dt * elem<double>{ r_tt, u_tt, v_tt, w_tt, p_tt } * (consts.num1 / consts.num12);
        get_elem(res.h, ix, iy, iz - Ns)
            = get_elem(buf.h, ix, iy, iz) + consts.dt * elem<double>{ r_t, u_t, v_t, w_t, p_t };

        asm volatile("// step end");
    }
}

void next(
    floors<double>* fio, /// 入出力
    tmp_walls<double>* tmp_wall, /// 壁バッファ
    local_buf<double>* tmp1, /// 一回で処理する分を入れるためのバッファ
    local_buf<double>* tmp2, /// 一回で処理する分を入れるためのバッファ
    int phase /// 0: 内部, 1: 皮
)
{
    TakeProfile<false, false> prof_all_(auxPerf[0]);

    const int tid = get_tid();
    const int gid = get_pid() * THREADS_PER_PE + tid;
    const int gws = THREADS_PER_PE * get_maxpid();

    elem2<double> e;

    // 領域に順にupdateを適用し、壁と床出力を受けとる
    for (int jx = MX - 1; jx >= (phase ? 0 : 1); jx--) {
        for (int jy = MY - 1; jy >= (phase ? 0 : 1); jy--) {
            for (int jz = MZ - 1; jz >= (phase ? 0 : 1); jz--) {
                // 皮フェーズの処理
                if (phase && jx != 0 && jy != 0 && jz != 0) {
                    continue;
                }

                // 床を読む
                int ix0 = jx * NX;
                int iy0 = jy * NY;
                int iz0 = jz * NZ;

                // start of update
                // split_cubo(tmp1->p, tmp1->h, fio->v, { 0, 0, 0 }, { ix0, iy0, iz0 }, { NX, NY, NZ });

                if (tid == 0 && get_pid() < NX * NY) {
                    // TakeProfile<false, true> prof_all_(auxPerf[5]);

                    asm volatile("// floor split start");

                    int x = get_pid() / NY;
                    int y = get_pid() % NY;

                    for (int z = 0; z < NZ; z++) {
                        const auto e = get_elem(fio->v, x + ix0, y + iy0, z + iz0);
                        get_elem(tmp1->p, x, y, z) = e.p;
                        get_elem(tmp1->h, x, y, z) = e.h;
                    }

                    asm volatile("// floor split end");
                }

                // flush_L2();
                // sync();
                // sync_L2();

                for (int it = 0; it < NT; it++) {
                    // x方向の壁を読み書き
                    // split_cubo(tmp1->p, tmp1->h, tmp_wall->x[jy][jz][it],
                    //     { NX, 0, 0 }, { 0, 0, 0 }, { 2 * Ns, NY + 2 * Ns, NZ + 2 * Ns });
                    // merge_cubo(tmp_wall->x[jy][jz][it], tmp1->p, tmp1->h,
                    //     { 0, 0, 0 }, { 0, 0, 0 }, { 2 * Ns, NY + 2 * Ns, NZ + 2 * Ns });
                    {
                        // TakeProfile<false, true> prof_all_(auxPerf[2]);
                        const int wx = 2 * Ns;
                        const int wy = 2 * Ns + NY;
                        const int wz = 2 * Ns + NZ;

                        if (tid == 0 && get_pid() < wy * wz) {
                            auto& tw = tmp_wall->x[jy][jz][it];
                            int y = get_pid() % wy;
                            int z = get_pid() / wy;

                            asm volatile("// xwall start");

                            auto p1 = tw[z][y];
                            auto p2 = tmp1->p[z][y];
                            auto p3 = tmp1->h[z][y];

                            for (int x = 0; x < Ns * 2; x++) {
                                const auto e = get_elem(tw, x + NX, y, z);
                                get_elem(tmp1->p, x, y, z) = e.p;
                                get_elem(tmp1->h, x, y, z) = e.h;
                                get_elem(tw, x, y, z) = e;
                            }
                            asm volatile("// xwall end");
                        }
                    }
                    // sync_L2();

                    // y方向の壁を読み書き
                    // split_cubo(tmp1->p, tmp1->h, tmp_wall->y[jz][jx][it],
                    //     { 0, NY, 0 }, { 0, 0, 0 }, { NX + 2 * Ns, 2 * Ns, NZ + 2 * Ns });
                    // merge_cubo(tmp_wall->y[jz][jx][it], tmp1->p, tmp1->h,
                    //     { 0, 0, 0 }, { 0, 0, 0 }, { NX + 2 * Ns, 2 * Ns, NZ + 2 * Ns });
                    {
                        // TakeProfile<false, true> prof_all_(auxPerf[3]);
                        const int wx = 2 * Ns + NX;
                        const int wy = 2 * Ns;
                        const int wz = 2 * Ns + NZ;

                        if (tid == 0 && get_pid() < wz * wx) {
                            auto& tw = tmp_wall->y[jz][jx][it];
                            int z = get_pid() % wz;
                            int x = get_pid() / wz;

                            for (int y = 0; y < Ns * 2; y++) {
                                const auto e = get_elem(tw, x, y + NY, z);
                                get_elem(tmp1->p, x, y, z) = e.p;
                                get_elem(tmp1->h, x, y, z) = e.h;
                                get_elem(tw, x, y, z) = e;
                            }
                        }
                    }
                    // sync_L2();

                    // z方向の壁を読み書き
                    // split_cubo(tmp1->p, tmp1->h, tmp_wall->z[jx][jy][it],
                    //     { 0, 0, NZ }, { 0, 0, 0 }, { NX + 2 * Ns, NY + 2 * Ns, 2 * Ns });
                    // merge_cubo(tmp_wall->z[jx][jy][it], tmp1->p, tmp1->h,
                    //     { 0, 0, 0 }, { 0, 0, 0 }, { NX + 2 * Ns, NY + 2 * Ns, 2 * Ns });
                    {
                        // TakeProfile<false, true> prof_all_(auxPerf[4]);
                        const int wx = 2 * Ns + NX;
                        const int wy = 2 * Ns + NY;
                        const int wz = 2 * Ns;
                        if (tid == 0 && get_pid() < wx * wy) {
                            auto& tw = tmp_wall->z[jx][jy][it];
                            int x = get_pid() / wy;
                            int y = get_pid() % wy;

                            for (int z = 0; z < Ns * 2; z++) {
                                const auto e = get_elem(tw, x, y, z + NZ);
                                get_elem(tmp1->p, x, y, z) = e.p;
                                get_elem(tmp1->h, x, y, z) = e.h;
                                get_elem(tw, x, y, z) = e;
                            }
                        }
                    }

                    sync();
                    flush_L2();

                    for (int cur = gid; cur < NX * NY * 4; cur += gws) {
                        int bz = gid / (NX * NY);
                        int tt = gid % (NX * NY);
                        int ix = (tt / 16) / (NX / 4) * 4 + (tt % 16 % 4);
                        int iy = (tt / 16) % (NX / 4) * 4 + (tt % 16 / 4);
                        sync_L2();
                        step(ix + Ns, iy + Ns, bz, NZ / 4, *tmp1, *tmp2);
                    }

                    flush_L2();
                    sync();
                    swap(tmp1, tmp2);
                }

                // merge_cubo(fio->v, tmp1->p, tmp1->h, { ix0, iy0, iz0 }, { 0, 0, 0 }, { NX, NY, NZ });
                if (tid == 0 && get_pid() < NX * NY) {
                    // TakeProfile<false, true> prof_all_(auxPerf[6]);

                    int x = get_pid() / NY;
                    int y = get_pid() % NY;

                    for (int z = 0; z < NZ; z++) {
                        auto& e = get_elem(fio->v, x + ix0, y + iy0, z + iz0);
                        e.p = get_elem(tmp1->p, x, y, z);
                        e.h = get_elem(tmp1->h, x, y, z);
                    }
                }

                sync_L2();
                // end of update
            }
        }
    }
}
