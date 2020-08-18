#include "update_sc.h"
#include "cache_profiler.h"
#include "config.h"
#include "pzcperf.h"

#include <mpi.h>

#include <assert.h>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class StopWatch {
public:
    StopWatch()
        : start(std::chrono::system_clock::now())
    {
    }

    double elapsed()
    {
        auto now = std::chrono::system_clock::now();
        return std::chrono::duration<double>(now - start).count();
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start;
};

static auto clExtSetPerThreadStackSize = (pfnPezyExtSetPerThreadStackSize)
    clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_set_per_thread_stack_size");

static auto clExtSetCacheWriteBuffer = (pfnPezyExtSetCacheWriteBuffer)
    clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_set_cache_writebuffer");

template <typename T, size_t X1, size_t Y1, size_t Z1, size_t X2, size_t Y2, size_t Z2>
void copy_rect(
    elem<T> (&dst)[X1][Y1][Z1],
    elem<T> (&src)[X2][Y2][Z2],
    const array<size_t, 3>& dest_offset,
    const array<size_t, 3>& src_offset,
    const array<size_t, 3>& region)
{
    for (int ix = 0; ix < region[0]; ix++) {
        for (int iy = 0; iy < region[1]; iy++) {
            for (int iz = 0; iz < region[2]; iz++) {
                dst[ix + dest_offset[0]][iy + dest_offset[1]][iz + dest_offset[2]] = src[ix + src_offset[0]][iy + src_offset[1]][iz + src_offset[2]];
            }
        }
    }
}

extern uint8_t _binary_kernel_pz_start;
extern uint8_t _binary_kernel_pz_end;

Updater::Updater(int device_id, int x_procs, int y_procs, int z_procs)
    : time_step_(0)
    , x_procs_(x_procs)
    , y_procs_(y_procs)
    , z_procs_(z_procs)
{
    assert((x_procs == 1 && y_procs == 1 && z_procs == 1) || (x_procs % 2 == 0 && y_procs % 2 == 0 && z_procs % 2 == 0));

    cl::Platform::get(&platform_);

    {
        vector<cl::Device> devices;
        platform_.getDevices(CL_DEVICE_TYPE_DEFAULT, &devices);
        assert(device_id < (int)devices.size());
        device_ = devices[device_id];
    }

    context_ = cl::Context(device_);
    queue_ = cl::CommandQueue(context_, device_);
    queue_trans_ = cl::CommandQueue(context_, device_);

    vector<uint8_t> kernel_source(&_binary_kernel_pz_start, &_binary_kernel_pz_end);

    try {
        program_ = cl::Program(context_, { device_ }, vector<vector<uint8_t>>{ kernel_source });
        program_.build({ device_ }, "--release");
    } catch (const cl::Error& e) {
        cout << e.what() << ": " << e.err() << endl;
        throw e;
    }

    clExtSetProfile(context_(), 0, true);
    clExtSetCacheWriteBuffer(context_(), 0, false);

    if (0) {
        uint64_t sps[8];
        cl::Buffer buf(context_, CL_MEM_READ_WRITE, sizeof(sps));
        cl::KernelFunctor<cl::Buffer&> k(program_, "get_stack_ptrs");

        // 2.5KB / thread
        assert(clExtSetPerThreadStackSize != nullptr);
        auto ret = clExtSetPerThreadStackSize(k.getKernel()(), 1024 * 5 / 2);
        printf("clExtSetPerThreadStackSize: %d\n", ret);
        fflush(stdout);

        k(cl::EnqueueArgs(queue_, 124 * 16 * 8), buf).wait();
        queue_.enqueueReadBuffer(buf, true, 0, sizeof(sps), sps);

        for (int i = 0; i < 8; i++) {
            printf("SP[%d] = %016lX\n", i, sps[i]);
        }
        fflush(stdout);

        exit(0);
    }

    if (0) {
        unique_ptr<local_buf<double>> p(new local_buf<double>);
        cl::Buffer buf(context_, CL_MEM_READ_WRITE, sizeof(*p));
        cl::KernelFunctor<cl::Buffer&> k(program_, "test1");

        // 2.5KB / thread
        assert(clExtSetPerThreadStackSize != nullptr);
        auto ret = clExtSetPerThreadStackSize(k.getKernel()(), 1024 * 5 / 2);
        printf("clExtSetPerThreadStackSize: %d\n", ret);
        fflush(stdout);

        k(cl::EnqueueArgs(queue_, 124 * 16 * 8), buf).wait();
        queue_.enqueueReadBuffer(buf, true, 0, sizeof(*p), &*p);

        for (int x = 0; x < NX + 2 * Ns; x++) {
            for (int y = 0; y < NY + 2 * Ns; y++) {
                for (int z = 0; z < NZ + 2 * Ns; z++) {
                    auto& pp = p->p[x][y][z];
                    auto& qq = p->h[x][y][z];
                    printf("P[%d][%d][%d] = {%f, %f, %f, %f, %f}\n", x, y, z, pp.r, pp.u, pp.v, pp.w, pp.p);
                    printf("P[%d][%d][%d] = {%f, %f, %f, %f, %f}\n", x, y, z, qq.r, qq.u, qq.v, qq.w, qq.p);
                }
            }
        }
        fflush(stdout);
        exit(0);
    }

    size_t total_device_mem_size = 0;
    floor_ = cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(floors<double>));
    total_device_mem_size += sizeof(floors<double>);

    tmp_wall_ = cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(tmp_walls<double>));
    total_device_mem_size += sizeof(tmp_walls<double>);

    local_buf_[0] = cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(local_buf<double>));
    local_buf_[1] = cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(local_buf<double>));
    total_device_mem_size += sizeof(local_buf<double>) * 2;

    // MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

    assert(x_procs_ * y_procs_ * z_procs_ == world_size_);

    if (world_size_ == 1) {
        x_pos_ = 0;
        y_pos_ = 0;
        z_pos_ = 0;
    } else {
        int x_cell = (my_rank_ / 8) / (z_procs_ / 2) / (y_procs_ / 2);
        int y_cell = (my_rank_ / 8) / (z_procs_ / 2) % (y_procs_ / 2);
        int z_cell = (my_rank_ / 8) % (z_procs_ / 2);

        x_pos_ = x_cell * 2 + (my_rank_ % 8 / 2 / 2);
        y_pos_ = y_cell * 2 + (my_rank_ % 8 / 2 % 2);
        z_pos_ = z_cell * 2 + (my_rank_ % 8 % 2);

        printf("[%03d]: (%d, %d, %d)\n", my_rank_, x_pos_, y_pos_, z_pos_);
    }

    if (my_rank_ == 0) {
        printf("Start calculating:\n");
        printf("  * NX=%d, NY=%d, NZ=%d, MX=%d, MY=%d, MZ=%d, NT=%d\n", NX, NY, NZ, MX, MY, MZ, NT);
        printf("  * Used device memory:   %.3fGB\n", (double)total_device_mem_size / (1 << 30));
        printf("  * Used host memory:     %.3fGB\n", (double)sizeof(Updater) / (1 << 30));
        printf("  * Size of global floor: %.3fGB\n", (double)sizeof(floors<double>) / (1 << 30));
        printf("  * Size of local floor:  %.3fMB\n", (double)sizeof(local_buf<double>) / (1 << 20));
        printf("  * Size of wall buffer:  %.3fGB\n", (double)sizeof(tmp_walls<double>) / (1 << 30));
        printf("  * FLOP per elem:        %d\n", FLOP_PER_ELEM);
        printf("  * FLOP per step:        %.3fG\n", (double)FLOP_PER_ELEM * LX * LY * LZ / 1.0e9);
        printf("  * FLOP per NT step:     %.3fG\n", (double)FLOP_PER_ELEM * LX * LY * LZ * NT / 1.0e9);
    }

    {
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

        printf("Size of diffs: %d, %d\n", sizeof(diffs), sizeof(diffs[5]));
    }
}

template <typename T, size_t AX, size_t AY, size_t AZ>
void read_cuboid(
    cl::CommandQueue& q,
    T (&arr)[AX][AY][AZ],
    cl::Buffer& buf,
    const array<size_t, 3>& buf_ofs,
    const array<size_t, 3>& region,
    size_t buf_ny, size_t buf_nz)
{
    q.enqueueReadBufferRect(buf, true,
        { buf_ofs[2] * sizeof(T), buf_ofs[1], buf_ofs[0] },
        { 0, 0, 0 },
        { region[2] * sizeof(T), region[1], region[0] },
        buf_nz * sizeof(T),
        buf_ny * buf_nz * sizeof(T),
        sizeof(arr[0][0]),
        sizeof(arr[0]),
        arr);
}

template <typename T, size_t AX, size_t AY, size_t AZ>
void write_cuboid(
    cl::CommandQueue& q,
    const T (&arr)[AX][AY][AZ],
    cl::Buffer& buf,
    const array<size_t, 3>& buf_ofs,
    const array<size_t, 3>& region,
    size_t buf_ny, size_t buf_nz)
{
    q.enqueueWriteBufferRect(buf, true,
        { buf_ofs[2] * sizeof(T), buf_ofs[1], buf_ofs[0] },
        { 0, 0, 0 },
        { region[2] * sizeof(T), region[1], region[0] },
        buf_nz * sizeof(T),
        buf_ny * buf_nz * sizeof(T),
        sizeof(arr[0][0]),
        sizeof(arr[0]),
        arr);
}

void Updater::transfer()
{
    vector<MPI_Request> reqs;

    const size_t buf_ny = MY * NY;
    const size_t buf_nz = MZ * NZ;

    const size_t thick = 2 * NT * Ns;

    // 面
    read_cuboid(queue_trans_, send_.x, floor_, { LX, 0, 0 }, { thick, LY, LZ }, buf_ny, buf_nz);
    read_cuboid(queue_trans_, send_.y, floor_, { 0, LY, 0 }, { LX, thick, LZ }, buf_ny, buf_nz);
    read_cuboid(queue_trans_, send_.z, floor_, { 0, 0, LZ }, { LX, LY, thick }, buf_ny, buf_nz);

    // 柱
    read_cuboid(queue_trans_, send_.xy, floor_, { LX, LY, 0 }, { thick, thick, LZ }, buf_ny, buf_nz);
    read_cuboid(queue_trans_, send_.yz, floor_, { 0, LY, LZ }, { LX, thick, thick }, buf_ny, buf_nz);
    read_cuboid(queue_trans_, send_.zx, floor_, { LX, 0, LZ }, { thick, LY, thick }, buf_ny, buf_nz);

    // 角
    read_cuboid(queue_trans_, send_.xyz, floor_, { LX, LY, LZ }, { thick, thick, thick }, buf_ny, buf_nz);

#define SENDRECV(buf, tag, dx, dy, dz)                                        \
    {                                                                         \
        MPI_Request req_send, req_recv;                                       \
        MPI_Isend(&send_.buf, sizeof(send_.buf) / sizeof(double), MPI_DOUBLE, \
            get_rank(x_pos_ + (dx), y_pos_ + (dy), z_pos_ + (dz)), (tag),     \
            MPI_COMM_WORLD, &req_send);                                       \
        MPI_Irecv(&recv_.buf, sizeof(recv_.buf) / sizeof(double), MPI_DOUBLE, \
            get_rank(x_pos_ - (dx), y_pos_ - (dy), z_pos_ - (dz)), (tag),     \
            MPI_COMM_WORLD, &req_recv);                                       \
        reqs.push_back(req_send);                                             \
        reqs.push_back(req_recv);                                             \
    }

    SENDRECV(x, 1100, 1, 0, 0);
    SENDRECV(y, 1010, 0, 1, 0);
    SENDRECV(z, 1001, 0, 0, 1);
    SENDRECV(xy, 1110, 1, 1, 0);
    SENDRECV(yz, 1011, 0, 1, 1);
    SENDRECV(zx, 1101, 1, 0, 1);
    SENDRECV(xyz, 1111, 1, 1, 1);

    vector<MPI_Status> stats(reqs.size());
    MPI_Waitall(reqs.size(), reqs.data(), stats.data());

    // 面
    write_cuboid(queue_trans_, recv_.y, floor_, { 0, 0, 0 }, { LX, thick, LZ }, buf_ny, buf_nz);
    write_cuboid(queue_trans_, recv_.z, floor_, { 0, 0, 0 }, { LX, LY, thick }, buf_ny, buf_nz);

    // 柱
    write_cuboid(queue_trans_, recv_.xy, floor_, { 0, 0, 0 }, { thick, thick, LZ }, buf_ny, buf_nz);
    write_cuboid(queue_trans_, recv_.yz, floor_, { 0, 0, 0 }, { LX, thick, thick }, buf_ny, buf_nz);
    write_cuboid(queue_trans_, recv_.zx, floor_, { 0, 0, 0 }, { thick, LY, thick }, buf_ny, buf_nz);

    // 角
    write_cuboid(queue_trans_, recv_.xyz, floor_, { 0, 0, 0 }, { thick, thick, thick }, buf_ny, buf_nz);
}

#define KERNEL_PROF
#define MAX_PERF 16384

void Updater::next()
{
    StopWatch sw_total;

    // カーネル呼び出し
    // 全体をNTステップ一気に更新する
    cl::KernelFunctor<
        cl::Buffer&, // floor
        cl::Buffer&, // tmp_wall
        cl::Buffer&, // local_buf
        cl::Buffer&, // local_buf
        int // phase
        >
        k(program_, "next");

    // 2.5KB/thread に設定
    assert(clExtSetPerThreadStackSize != nullptr);
    clExtSetPerThreadStackSize(k.getKernel()(), 1024 * 5 / 2);

    if (my_rank_ == 0) {
        printf("Run kernel phase 0\n");
        fflush(stdout);
    }

#ifdef KERNEL_PROF
    Profiles prof_start(context_);
#endif

    StopWatch sw_kernel;

    // 皮部分以外の計算を開始
    auto ev = k(cl::EnqueueArgs(queue_, 124 * 16 * 8),
        floor_, tmp_wall_, local_buf_[0], local_buf_[1], 0);

    // 皮部分以外を実行中に転送を終わらせる
    {
        StopWatch sw;
        transfer();
        if (my_rank_ == 0) {
            printf("Transfer: %.3fms\n", sw.elapsed() * 1e3);
        }
    }

    ev.wait();

    if (my_rank_ == 0) {
        printf("Kernel phase 0: %.3fms\n", sw_kernel.elapsed() * 1e3);
    }

#ifdef KERNEL_PROF
    {
        PZCPerformance perf[16];
        cl::KernelFunctor<cl::Buffer&, int> getperf(program_, "GetPerformance");
        cl::Buffer buf[16];
        for (int i = 0; i < 16; i++) {
            buf[i] = cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(PZCPerformance));
        }

        for (int i = 0; i < 16; i++) {
            getperf(cl::EnqueueArgs(queue_, 124 * 16 * 8), buf[i], MAX_PERF + i).wait();
            queue_.enqueueReadBuffer(buf[i], true, 0, sizeof(perf[i]), &perf[i]);
        }

        auto print_perf = [](const PZCPerformance& perf) {
            printf("perf = %.3fM, stall = %.3fM, istall = %.3fM, wait = %.3fM\n",
                perf.perf / 1e6, perf.stall / 1e6, perf.istall / 1e6, perf.wait / 1e6);
        };

        printf("[all]: ");
        print_perf(perf[0]);

        printf("[stp]: ");
        print_perf(perf[1]);

        perf[1] = perf[0] - perf[1];
        printf("[oth]: ");
        print_perf(perf[1]);

        for (int i = 2; i < 16; i++) {
            if (perf[i].perf != 0) {
                printf("[00%d]: ", i);
                print_perf(perf[i]);
            }
        }
    }
#endif

    // 皮計算
    {
        StopWatch sw;
        k(cl::EnqueueArgs(queue_, 124 * 16 * 8),
            floor_, tmp_wall_, local_buf_[0], local_buf_[1], 1)
            .wait();
        if (my_rank_ == 0) {
            printf("Kernel phase 1: %.3fms\n", sw.elapsed() * 1e3);
        }
    }

#ifdef KERNEL_PROF
    if (my_rank_ == 0) {
        Profiles prof_end(context_);
        prof_end -= prof_start;

        prof_end.print();

        double flop_per_elem = 2829;
        double flop_per_process = (double)NT * LX * LY * LZ * flop_per_elem;
        double flop = flop_per_process * world_size_;
        double sec = prof_end.pe_stats.elapse_ns / 1e9;
        double flops = flop / sec;
        double eff = prof_end.pe_stats.efficiency;

        printf("%.3f GFlops = %.3f GFlop / %.3f msec [efficiency = %.3f%%]\n",
            flops / 1e9, flop / 1e9, sec * 1e3, eff);
        fflush(stdout);
    }

    if (my_rank_ == 0) {
        printf("Total: %.3fms\n", sw_total.elapsed() * 1e3);
    }
#endif

#if 0

  // tmp_wall_xの初期化
  elem<double> init = {1, 0, 0, 0, 0};

  for (int jz = 0; jz < MZ+2; jz++) {
    for (int jy = 0; jy < MY+2; jy++) {
      for(int it = 0; it < NT; it++) {
        for(int ix = 0; ix < 2*Ns; ix++){
          for(int iy = 0; iy < NY+2*Ns; iy++){
            for(int iz = 0; iz < NZ+2*Ns; iz++){
              p_tmp_wall_x[jy][jz][it][ix][iy][iz] = init;
              h_tmp_wall_x[jy][jz][it][ix][iy][iz] = init;
            }
          }
        }
      }
    }
  }

  // 領域に順にupdateを適用し、壁と床出力を受けとる
  for(int jx = MX+1;jx >= 0; jx--) {
    // tmp_wall_yの初期化
    for (int jz = 0; jz < MZ+2; jz++) {
      for(int it = 0; it < NT; it++){
        for(int ix = 0; ix < NX+2*Ns; ix++){
          for(int iy = 0; iy < 2*Ns; iy++){
            for(int iz = 0; iz < NZ+2*Ns; iz++){
              p_tmp_wall_y[jz][it][ix][iy][iz] = init;
              h_tmp_wall_y[jz][it][ix][iy][iz] = init;
            }
          }
        }
      }
    }

    for(int jy = MY+1;jy >= 0; jy--) {
      // z壁の初期化
      for(int it = 0; it < NT; it++){
        for(int ix = 0; ix < NX+2*Ns; ix++){
          for(int iy = 0; iy < NY+2*Ns; iy++){
            for(int iz = 0; iz < 2*Ns; iz++){
              p_wall_z[it][ix][iy][iz] = init;
              h_wall_z[it][ix][iy][iz] = init;
            }
          }
        }
      }

      for(int jz = MZ+1;jz >= 0; jz--) {
        x_origin = jx*NX;
        y_origin = jy*NY;
        z_origin = jz*NZ;

        for(int it = 0; it < NT; it++) {
          // x壁を渡す
          for(int ix = 0; ix < 2*Ns; ix++){
            for(int iy = 0; iy < NY+2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p_wall_x[it][ix][iy][iz] = p_tmp_wall_x[jy][jz][it][ix][iy][iz];
                h_wall_x[it][ix][iy][iz] = h_tmp_wall_x[jy][jz][it][ix][iy][iz];
              }
            }
          }

          // y壁を渡す
          for(int ix = 0; ix < NX+2*Ns; ix++){
            for(int iy = 0; iy < 2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p_wall_y[it][ix][iy][iz] = p_tmp_wall_y[jz][it][ix][iy][iz];
                h_wall_y[it][ix][iy][iz] = h_tmp_wall_y[jz][it][ix][iy][iz];
              }
            }
          }
        }

        update(s);

        for(int it = 0; it < NT; it++) {
          // x壁の受けとり
          for(int ix = 0; ix < 2*Ns; ix++){
            for(int iy = 0; iy < NY+2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p_tmp_wall_x[jy][jz][it][ix][iy][iz] = p_wall_x[it][ix][iy][iz];
                h_tmp_wall_x[jy][jz][it][ix][iy][iz] = h_wall_x[it][ix][iy][iz];
              }
            }
          }

          // y壁の受けとり
          for(int ix = 0; ix < NX+2*Ns; ix++){
            for(int iy = 0; iy < 2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p_tmp_wall_y[jz][it][ix][iy][iz] = p_wall_y[it][ix][iy][iz];
                h_tmp_wall_y[jz][it][ix][iy][iz] = h_wall_y[it][ix][iy][iz];
              }
            }
          }
        }

        // 床の受けとり
        for(int ix = 0; ix < NX; ix++) {
          for(int iy = 0; iy < NY; iy++) {
            for(int iz = 0; iz < NZ; iz++) {
              p_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p_res[ix][iy][iz];
              h_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = h_res[ix][iy][iz];
            }
          }
        }
      }
    }
  }

  // 床を書き込む
  for(int ix = 0; ix < LX; ix++) {
    int x = (ix - NX + LX)%LX + 2*NX;
    /* int x = ix + 2*NX; */
    for(int iy = 0; iy < LY; iy++) {
      int y = (iy - NY + LY)%LY + 2*NY;
      /* int y = iy + 2*NY; */
      for(int iz = 0; iz < LZ; iz++) {
        int z = (iz - NZ + LZ)%LZ + 2*NZ;
        /* int z = iz + 2*NZ; */
        s.p[x][y][z] = p_tmp[ix][iy][iz];
        s.h[x][y][z] = h_tmp[ix][iy][iz];
      }
    }
  }

#endif

    // NTだけタイムステップを進める
    time_step_ += NT;
}
