#include "config.h"
#include "update_sc.h"

#include "cmdline.h"

#include <mpi.h>

#include <math.h>
#include <stdio.h>
#include <cassert>
#include <memory>
#include <chrono>

using namespace std;

double x0 = LX*h/2.0;

void write_data(const navi &n, const floors<double> &s) {
  double t;
  char fn[256];
  char fn_xy[256];
  char fn_xz[256];

  t = n.time_step * dt;
  printf("it = %d: t = %f\n", n.time_step, t);

  sprintf(fn, "data/%f.dat", t);
  sprintf(fn_xy, "data/%f-xy.dat", t);
  sprintf(fn_xz, "data/%f-xz.dat", t);
  /* sprintf(fn, "data/%f_%d.dat", t, mpi_my_rank); */
  /* sprintf(fn_xy, "data/%f_%d-xy.dat", t, mpi_my_rank); */
  /* sprintf(fn_xz, "data/%f_%d-xz.dat", t, mpi_my_rank); */
  FILE *fp = fopen(fn, "w");
  FILE *fp_xy = fopen(fn_xy, "w");
  FILE *fp_xz = fopen(fn_xz, "w");
  for(int ix = n.lower_x; ix < n.upper_x; ++ix) {
    double x = (ix + n.offset_x)*h;
    for(int iy = n.lower_y; iy < n.upper_y; ++iy) {
      double y = (iy + n.offset_y)*h;
      for(int iz = n.lower_z; iz < n.upper_z; ++iz) {
        double z = (iz + n.offset_z)*h;

// FIXME
/*
        fprintf(fp, "%.5f %.5f %.5f %f %f %f %f %f %f %f %f %f %f\n", x, y, z,
                s.p[ix][iy][iz].r,
                s.p[ix][iy][iz].u,
                s.p[ix][iy][iz].v,
                s.p[ix][iy][iz].w,
                s.p[ix][iy][iz].p,
                s.h[ix][iy][iz].r,
                s.h[ix][iy][iz].u,
                s.h[ix][iy][iz].v,
                s.h[ix][iy][iz].w,
                s.h[ix][iy][iz].p);

        if (z == x0) {
          fprintf(fp_xy, "%.5f %.5f %.5f %f %f %f %f %f %f %f %f %f %f\n", x, y, z,
                  s.p[ix][iy][iz].r,
                  s.p[ix][iy][iz].u,
                  s.p[ix][iy][iz].v,
                  s.p[ix][iy][iz].w,
                  s.p[ix][iy][iz].p,
                  s.h[ix][iy][iz].r,
                  s.h[ix][iy][iz].u,
                  s.h[ix][iy][iz].v,
                  s.h[ix][iy][iz].w,
                  s.h[ix][iy][iz].p);
        }

        if (y == x0) {
          fprintf(fp_xz, "%.5f %.5f %.5f %f %f %f %f %f %f %f %f %f %f\n", x, y, z,
                  s.p[ix][iy][iz].r,
                  s.p[ix][iy][iz].u,
                  s.p[ix][iy][iz].v,
                  s.p[ix][iy][iz].w,
                  s.p[ix][iy][iz].p,
                  s.h[ix][iy][iz].r,
                  s.h[ix][iy][iz].u,
                  s.h[ix][iy][iz].v,
                  s.h[ix][iy][iz].w,
                  s.h[ix][iy][iz].p);
        }
*/
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp_xy, "\n");
    fprintf(fp_xz, "\n");
  }
  fclose(fp);
  fclose(fp_xy);
  fclose(fp_xz);
  printf("write: data/%f.dat\n", t);

  /* printf("write: data/%f_%d.dat\n", t, mpi_my_rank); */
}

void check_conservation(const navi &n, const floors<double> &s) {
// FIXME
/*
  double mp = 0.0;
  double mh = 0.0;
  double ekp = 0.0;
  double ekh = 0.0;
  double eip = 0.0;
  double eih = 0.0;
  double ep = 0.0;
  double eh = 0.0;
  for(int ix = n.lower_x; ix < n.upper_x; ++ix) {
    for(int iy = n.lower_y; iy < n.upper_y; ++iy) {
      for(int iz = n.lower_z; iz < n.upper_z; ++iz) {
        double rp0 = s.p[ix][iy][iz].r;
        double up0 = s.p[ix][iy][iz].u;
        double vp0 = s.p[ix][iy][iz].v;
        double wp0 = s.p[ix][iy][iz].w;
        double pp0 = s.p[ix][iy][iz].p;
        double rh0 = s.h[ix][iy][iz].r;
        double uh0 = s.h[ix][iy][iz].u;
        double vh0 = s.h[ix][iy][iz].v;
        double wh0 = s.h[ix][iy][iz].w;
        double ph0 = s.h[ix][iy][iz].p;

        double vol = h*h*h;

        mp += rp0*vol;
        mh += rh0*vol;
        ekp += (rp0*(up0*up0 + vp0*vp0 + wp0*wp0)/2.0)*vol;
        ekh += (rh0*(uh0*uh0 + vh0*vh0 + wh0*wh0)/2.0)*vol;
        eip += vol*pp0/(gm-1);
        eih += vol*ph0/(gm-1);
        ep += (rp0*(up0*up0 + vp0*vp0 + wp0*wp0)/2.0 + pp0/(gm-1))*vol;
        eh += (rh0*(uh0*uh0 + vh0*vh0 + wh0*wh0)/2.0 + ph0/(gm-1))*vol;
      }
    }
  }
  printf("P: m = %f, e = %f + %f = %f\n", mp, ekp, eip, ep);
  printf("H: m = %f, e = %f + %f = %f\n", mh, ekh, eih, eh);
*/
}

/*
void main_cpu() {
  navi n;
  state st;

  init(n, st);

  printf("NX = %d\n", NX);
  printf("MX = %d\n", MX);
  printf("NT = %d\n", NT);
  printf("h = %f\n", h);
  printf("dt = %f\n", dt);
  printf("x0 = %f\n", x0);

  while(n.time_step < T_MAX) {
    fflush(stdout);

    write_data(n, st);
    check_conservation(n, st);
    next(n, st);

    if (n.time_step % FILTER_INTERVAL == 0) {
      lowpass_filter(n, st);
      printf("filtered\n");
    }
  }
  write_data(n, st);
  check_conservation(n, st);
}
*/

void main_sc(int px, int py, int pz) {
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int device_id = rank % 8;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  assert(size == px * py * pz);

  navi n;
  auto upd = unique_ptr<Updater>(new Updater(device_id, px, py, pz));
  // auto st = unique_ptr<state>(new state());
  // init(n, *st);

  // printf("NX = %d\n", NX);
  // printf("MX = %d\n", MX);
  // printf("NT = %d\n", NT);
  // printf("h = %f\n", h);
  // printf("dt = %f\n", dt);
  // printf("x0 = %f\n", x0);

  auto start = chrono::system_clock::now();

  while(upd->time_step() < T_MAX) {
    auto prev_step = upd->time_step();
    auto cur_start = chrono::system_clock::now();

    // write_data(n);
    // check_conservation(n);
    upd->next();

    // if (n.time_step % FILTER_INTERVAL == 0) {
    //   lowpass_filter(n);
    //   printf("filtered\n");
    // }

    auto steps = upd->time_step() - prev_step;
    chrono::duration<double> elapsed = chrono::system_clock::now() - cur_start;
    chrono::duration<double> total = chrono::system_clock::now() - start;

    if (rank == 0) {
      uint64_t flop_per_step = (uint64_t)FLOP_PER_ELEM * LX * LY * LZ * size;

      printf("[%ld]: Current = %.3fGFlops, Average = %.3fGFlops (%.3fGFlops / Process)\n",
             upd->time_step(),
             (double)flop_per_step * steps / elapsed.count() / 1e9,
             (double)flop_per_step * upd->time_step() / total.count() / 1e9,
             (double)flop_per_step * upd->time_step() / total.count() / 1e9 / size);
      fflush(stdout);
    }
  }

  // write_data(n, *st);
  // check_conservation(n, *st);
}

int main(int argc, char **argv)
{
  cmdline::parser a;
  a.add<int>("px", '\0', "Process number of x-direction", true);
  a.add<int>("py", '\0', "Process number of y-direction", true);
  a.add<int>("pz", '\0', "Process number of z-direction", true);

  a.parse_check(argc, argv);

  auto px = a.get<int>("px");
  auto py = a.get<int>("py");
  auto pz = a.get<int>("pz");

  assert(px * py * pz % 8 == 0 || (px == 1 && py == 1 && pz == 1));

  MPI_Init(&argc, &argv);

  main_sc(px, py, pz);

  MPI_Finalize();
  return 0;
}
