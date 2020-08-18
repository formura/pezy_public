#pragma once

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#include <CL/cl2.hpp>

#include "config.h"

// typedef elem<double> buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];

template <typename T>
struct transfer_buf {
  elem2<T> x[2 * NT * Ns][LY][LZ];
  elem2<T> y[LX][2 * NT * Ns][LZ];
  elem2<T> z[LX][LY][2 * NT * Ns];

  elem2<T> xy[2 * NT * Ns][2 * NT * Ns][LZ];
  elem2<T> yz[LX][2 * NT * Ns][2 * NT * Ns];
  elem2<T> zx[2 * NT * Ns][LY][2 * NT * Ns];

  elem2<T> xyz[2 * NT * Ns][2 * NT * Ns][2 * NT * Ns];
};

class Updater {
public:
  explicit Updater(int device_id, int x_procs, int y_procs, int z_procs);

  void next();

  uint64_t time_step() const {
    return time_step_;
  }

private:
  void transfer();

  // void update();

  int get_rank(int x, int y, int z) const {
    x = (x + x_procs_) % x_procs_;
    y = (y + y_procs_) % y_procs_;
    z = (z + z_procs_) % z_procs_;

    int node = (x / 2) * (y_procs_ / 2) * (z_procs_ / 2) +
           (y / 2) * (z_procs_ / 2) + (z / 2);

    return node * 8 + (x % 2 * 4) + (y % 2 * 2) + z % 2;
  }

  uint64_t time_step_;

  // OpenCL stuffs
  cl::Platform platform_;
  cl::Device device_;
  cl::Context context_;
  cl::CommandQueue queue_, queue_trans_;
  cl::Program program_;

  cl::Buffer floor_;
  cl::Buffer tmp_wall_;
  cl::Buffer local_buf_[2];

  // MPI stuffs
  int my_rank_, world_size_;
  int x_procs_, y_procs_, z_procs_;
  int x_pos_, y_pos_, z_pos_;

  // 送受信バッファ
  transfer_buf<double> send_, recv_;
};
