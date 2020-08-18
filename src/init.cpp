#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <memory>
#include "config.h"

using namespace std;

// 初期条件の作成

#if defined(VORTEX)
void setup_vortex() {
  printf("Setup vortex\n");
  double lx = h*LX;
  double ly = h*LY;
  for(int ix = 0; ix < LX+2*Ns; ix++) {
    double x = ((ix-Ns+LX)%LX)*h;
    for(int iy = 0; iy < LY+2*Ns; iy++) {
      double y = ((iy-Ns+LY)%LY)*h;
      for(int iz = 0; iz < LZ+2*Ns; iz++) {
        rc[ix][iy][iz] = 1.0;
        uc[ix][iy][iz] = -sin(4*M_PI*y/ly);
        vc[ix][iy][iz] = sin(4*M_PI*x/lx);
        wc[ix][iy][iz] = 0.0;
        pc[ix][iy][iz] = 1.0;
      }
    }
  }
}
#else
struct tmp_buf {
  double rc[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns];
  double uc[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns];
  double vc[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns];
  double wc[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns];
  double pc[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns];
  double vel1[LX][LY][LZ];
  double vel2[LX][LY][LZ];
  double vel3[LX][LY][LZ];
};

double energy_spectra(double k) {
  return A*pow(k,4.0)*exp(-B*k*k);
}

// 一様等方乱流場の生成
// Ref: 5.2 in https://doi.org/10.1016/j.jcp.2012.10.005
void setup() {
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  fftw_complex *vx1, *vx2, *vx3;
  fftw_complex *vk1, *vk2, *vk3;
  vx1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);
  vx2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);
  vx3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);
  vk1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);
  vk2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);
  vk3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * LX * LY * LZ);

  // -1 < v < 1 の一様乱数を生成しvx1,vx2,vx3を初期化する
  for(int ix = 0; ix < LX; ix++) {
    for(int iy = 0; iy < LY; iy++) {
      for(int iz = 0; iz < LZ; iz++) {
        int idx = ix + iy*LX + iz*LX*LY;
        // real part
        vx1[idx][0] = 2.0*(gsl_rng_uniform_pos(r)-0.5);
        vx2[idx][0] = 2.0*(gsl_rng_uniform_pos(r)-0.5);
        vx3[idx][0] = 2.0*(gsl_rng_uniform_pos(r)-0.5);
        // imaginary part
        vx1[idx][1] = 0.0;
        vx2[idx][1] = 0.0;
        vx3[idx][1] = 0.0;
      }
    }
  }

  // フーリエ変換 (x->k)
  fftw_plan pf1, pf2, pf3;
  pf1 = fftw_plan_dft_3d(LX, LY, LZ, vx1, vk1, FFTW_FORWARD, FFTW_ESTIMATE);
  pf2 = fftw_plan_dft_3d(LX, LY, LZ, vx2, vk2, FFTW_FORWARD, FFTW_ESTIMATE);
  pf3 = fftw_plan_dft_3d(LX, LY, LZ, vx3, vk3, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(pf1);
  fftw_execute(pf2);
  fftw_execute(pf3);

  double lx = LX*h;
  double ly = LY*h;
  double lz = LZ*h;
  double v1, v2, v3;
  for(int jx = 0; jx < LX; jx++) {
    double kx = 2*M_PI*jx/lx;
    for(int jy = 0; jy < LY; jy++) {
      double ky = 2*M_PI*jy/ly;
      for(int jz = 0; jz < LZ; jz++) {
        double kz = 2*M_PI*jz/lz;
        int idx = jx + jy*LX + jz*LX*LY;
        double k = sqrt(kx*kx + ky*ky + kz*kz);
        if (1 <= k && k <= 8) {
          v1 = vk1[idx][0];
          v2 = vk2[idx][0];
          v3 = vk3[idx][0];

          // 非圧縮条件を満たすように変換
          double kv = kx*v1 + ky*v2 + kz*v3;
          v1 = v1 - kx*kv/(k*k);
          v2 = v2 - ky*kv/(k*k);
          v3 = v3 - kz*kv/(k*k);

          // エネルギー分布に従うように変換
          double v = sqrt(v1*v1 + v2*v2 + v3*v3);
          vk1[idx][0] = v1*sqrt(energy_spectra(k)/(4*M_PI*k*k))/v;
          vk2[idx][0] = v2*sqrt(energy_spectra(k)/(4*M_PI*k*k))/v;
          vk3[idx][0] = v3*sqrt(energy_spectra(k)/(4*M_PI*k*k))/v;
        } else {
          vk1[idx][0] = 0.0;
          vk2[idx][0] = 0.0;
          vk3[idx][0] = 0.0;
        }
      }
    }
  }


  // フーリエ逆変換 (k->x)
  // 正規化されていないので注意
  fftw_plan pb1, pb2, pb3;
  pb1 = fftw_plan_dft_3d(LX, LY, LZ, vk1, vx1, FFTW_BACKWARD, FFTW_ESTIMATE);
  pb2 = fftw_plan_dft_3d(LX, LY, LZ, vk2, vx2, FFTW_BACKWARD, FFTW_ESTIMATE);
  pb3 = fftw_plan_dft_3d(LX, LY, LZ, vk3, vx3, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(pb1);
  fftw_execute(pb2);
  fftw_execute(pb3);

  // 結果を出力
  for(int ix = 0; ix < LX; ix++) {
    for(int iy = 0; iy < LY; iy++) {
      for(int iz = 0; iz < LZ; iz++) {
        int idx = ix + iy*LX + iz*LX*LY;
        vel1[ix][iy][iz] = vx1[idx][0]/(LX*LY*LZ);
        vel2[ix][iy][iz] = vx2[idx][0]/(LX*LY*LZ);
        vel3[ix][iy][iz] = vx3[idx][0]/(LX*LY*LZ);
      }
    }
  }

  fftw_destroy_plan(pf1);
  fftw_destroy_plan(pf2);
  fftw_destroy_plan(pf3);
  fftw_destroy_plan(pb1);
  fftw_destroy_plan(pb2);
  fftw_destroy_plan(pb3);
  fftw_free(vx1);
  fftw_free(vx2);
  fftw_free(vx3);
  fftw_free(vk1);
  fftw_free(vk2);
  fftw_free(vk3);
  gsl_rng_free(r);
}

void setup_turbulence(tmp_buf &tmp) {
  printf("Setup turbulence\n");
  // vel1, vel2, vel3 の初期化
  setup();

  for(int ix = 0; ix < LX+2*Ns; ix++) {
    int x = (ix-Ns+LX)%LX;
    for(int iy = 0; iy < LY+2*Ns; iy++) {
      int y = (iy-Ns+LY)%LY;
      for(int iz = 0; iz < LZ+2*Ns; iz++) {
        int z = (iz-Ns+LZ)%LZ;
        tmp.rc[ix][iy][iz] = 1.0;
        tmp.uc[ix][iy][iz] = vel1[x][y][z];
        tmp.vc[ix][iy][iz] = vel2[x][y][z];
        tmp.wc[ix][iy][iz] = vel3[x][y][z];
        tmp.pc[ix][iy][iz] = 1.0;
      }
    }
  }
}
#endif

double d_x(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j+0][k+0] -8 * q[i-1][j+0][k+0] +8 * q[i+1][j+0][k+0] -1 * q[i+2][j+0][k+0])/(12*h);}
double d_y(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j-2][k+0] -8 * q[i+0][j-1][k+0] +8 * q[i+0][j+1][k+0] -1 * q[i+0][j+2][k+0])/(12*h);}
double d_z(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j+0][k-2] -8 * q[i+0][j+0][k-1] +8 * q[i+0][j+0][k+1] -1 * q[i+0][j+0][k+2])/(12*h);}
double d_xx(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k+0] +16 * q[i-1][j+0][k+0] -30 * q[i+0][j+0][k+0] +16 * q[i+1][j+0][k+0] -1 * q[i+2][j+0][k+0])/(12*h*h);}
double d_xy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j-2][k+0] -8 * q[i-2][j-1][k+0] +8 * q[i-2][j+1][k+0] -1 * q[i-2][j+2][k+0] -8 * q[i-1][j-2][k+0] +64 * q[i-1][j-1][k+0] -64 * q[i-1][j+1][k+0] +8 * q[i-1][j+2][k+0] +8 * q[i+1][j-2][k+0] -64 * q[i+1][j-1][k+0] +64 * q[i+1][j+1][k+0] -8 * q[i+1][j+2][k+0] -1 * q[i+2][j-2][k+0] +8 * q[i+2][j-1][k+0] -8 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h);}
double d_yy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k+0] +16 * q[i+0][j-1][k+0] -30 * q[i+0][j+0][k+0] +16 * q[i+0][j+1][k+0] -1 * q[i+0][j+2][k+0])/(12*h*h);}
double d_xz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j+0][k-2] -8 * q[i-2][j+0][k-1] +8 * q[i-2][j+0][k+1] -1 * q[i-2][j+0][k+2] -8 * q[i-1][j+0][k-2] +64 * q[i-1][j+0][k-1] -64 * q[i-1][j+0][k+1] +8 * q[i-1][j+0][k+2] +8 * q[i+1][j+0][k-2] -64 * q[i+1][j+0][k-1] +64 * q[i+1][j+0][k+1] -8 * q[i+1][j+0][k+2] -1 * q[i+2][j+0][k-2] +8 * q[i+2][j+0][k-1] -8 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h);}
double d_yz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j-2][k-2] -8 * q[i+0][j-2][k-1] +8 * q[i+0][j-2][k+1] -1 * q[i+0][j-2][k+2] -8 * q[i+0][j-1][k-2] +64 * q[i+0][j-1][k-1] -64 * q[i+0][j-1][k+1] +8 * q[i+0][j-1][k+2] +8 * q[i+0][j+1][k-2] -64 * q[i+0][j+1][k-1] +64 * q[i+0][j+1][k+1] -8 * q[i+0][j+1][k+2] -1 * q[i+0][j+2][k-2] +8 * q[i+0][j+2][k-1] -8 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h);}
double d_zz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j+0][k-2] +16 * q[i+0][j+0][k-1] -30 * q[i+0][j+0][k+0] +16 * q[i+0][j+0][k+1] -1 * q[i+0][j+0][k+2])/(12*h*h);}
double d_xxx(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k+0] +2 * q[i-1][j+0][k+0] -2 * q[i+1][j+0][k+0] +1 * q[i+2][j+0][k+0])/(2*h*h*h);}
double d_xxy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +8 * q[i-2][j-1][k+0] -8 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +16 * q[i-1][j-2][k+0] -128 * q[i-1][j-1][k+0] +128 * q[i-1][j+1][k+0] -16 * q[i-1][j+2][k+0] -30 * q[i+0][j-2][k+0] +240 * q[i+0][j-1][k+0] -240 * q[i+0][j+1][k+0] +30 * q[i+0][j+2][k+0] +16 * q[i+1][j-2][k+0] -128 * q[i+1][j-1][k+0] +128 * q[i+1][j+1][k+0] -16 * q[i+1][j+2][k+0] -1 * q[i+2][j-2][k+0] +8 * q[i+2][j-1][k+0] -8 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h);}
double d_xyy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +16 * q[i-2][j-1][k+0] -30 * q[i-2][j+0][k+0] +16 * q[i-2][j+1][k+0] -1 * q[i-2][j+2][k+0] +8 * q[i-1][j-2][k+0] -128 * q[i-1][j-1][k+0] +240 * q[i-1][j+0][k+0] -128 * q[i-1][j+1][k+0] +8 * q[i-1][j+2][k+0] -8 * q[i+1][j-2][k+0] +128 * q[i+1][j-1][k+0] -240 * q[i+1][j+0][k+0] +128 * q[i+1][j+1][k+0] -8 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -16 * q[i+2][j-1][k+0] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h);}
double d_yyy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k+0] +2 * q[i+0][j-1][k+0] -2 * q[i+0][j+1][k+0] +1 * q[i+0][j+2][k+0])/(2*h*h*h);}
double d_xxz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +8 * q[i-2][j+0][k-1] -8 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +16 * q[i-1][j+0][k-2] -128 * q[i-1][j+0][k-1] +128 * q[i-1][j+0][k+1] -16 * q[i-1][j+0][k+2] -30 * q[i+0][j+0][k-2] +240 * q[i+0][j+0][k-1] -240 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] +16 * q[i+1][j+0][k-2] -128 * q[i+1][j+0][k-1] +128 * q[i+1][j+0][k+1] -16 * q[i+1][j+0][k+2] -1 * q[i+2][j+0][k-2] +8 * q[i+2][j+0][k-1] -8 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h);}
double d_xyz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j-2][k-2] -8 * q[i-2][j-2][k-1] +8 * q[i-2][j-2][k+1] -1 * q[i-2][j-2][k+2] -8 * q[i-2][j-1][k-2] +64 * q[i-2][j-1][k-1] -64 * q[i-2][j-1][k+1] +8 * q[i-2][j-1][k+2] +8 * q[i-2][j+1][k-2] -64 * q[i-2][j+1][k-1] +64 * q[i-2][j+1][k+1] -8 * q[i-2][j+1][k+2] -1 * q[i-2][j+2][k-2] +8 * q[i-2][j+2][k-1] -8 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] -8 * q[i-1][j-2][k-2] +64 * q[i-1][j-2][k-1] -64 * q[i-1][j-2][k+1] +8 * q[i-1][j-2][k+2] +64 * q[i-1][j-1][k-2] -512 * q[i-1][j-1][k-1] +512 * q[i-1][j-1][k+1] -64 * q[i-1][j-1][k+2] -64 * q[i-1][j+1][k-2] +512 * q[i-1][j+1][k-1] -512 * q[i-1][j+1][k+1] +64 * q[i-1][j+1][k+2] +8 * q[i-1][j+2][k-2] -64 * q[i-1][j+2][k-1] +64 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] +8 * q[i+1][j-2][k-2] -64 * q[i+1][j-2][k-1] +64 * q[i+1][j-2][k+1] -8 * q[i+1][j-2][k+2] -64 * q[i+1][j-1][k-2] +512 * q[i+1][j-1][k-1] -512 * q[i+1][j-1][k+1] +64 * q[i+1][j-1][k+2] +64 * q[i+1][j+1][k-2] -512 * q[i+1][j+1][k-1] +512 * q[i+1][j+1][k+1] -64 * q[i+1][j+1][k+2] -8 * q[i+1][j+2][k-2] +64 * q[i+1][j+2][k-1] -64 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] -1 * q[i+2][j-2][k-2] +8 * q[i+2][j-2][k-1] -8 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] +8 * q[i+2][j-1][k-2] -64 * q[i+2][j-1][k-1] +64 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] -8 * q[i+2][j+1][k-2] +64 * q[i+2][j+1][k-1] -64 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h);}
double d_yyz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +8 * q[i+0][j-2][k-1] -8 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +16 * q[i+0][j-1][k-2] -128 * q[i+0][j-1][k-1] +128 * q[i+0][j-1][k+1] -16 * q[i+0][j-1][k+2] -30 * q[i+0][j+0][k-2] +240 * q[i+0][j+0][k-1] -240 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] +16 * q[i+0][j+1][k-2] -128 * q[i+0][j+1][k-1] +128 * q[i+0][j+1][k+1] -16 * q[i+0][j+1][k+2] -1 * q[i+0][j+2][k-2] +8 * q[i+0][j+2][k-1] -8 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h);}
double d_xzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +16 * q[i-2][j+0][k-1] -30 * q[i-2][j+0][k+0] +16 * q[i-2][j+0][k+1] -1 * q[i-2][j+0][k+2] +8 * q[i-1][j+0][k-2] -128 * q[i-1][j+0][k-1] +240 * q[i-1][j+0][k+0] -128 * q[i-1][j+0][k+1] +8 * q[i-1][j+0][k+2] -8 * q[i+1][j+0][k-2] +128 * q[i+1][j+0][k-1] -240 * q[i+1][j+0][k+0] +128 * q[i+1][j+0][k+1] -8 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -16 * q[i+2][j+0][k-1] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h);}
double d_yzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +16 * q[i+0][j-2][k-1] -30 * q[i+0][j-2][k+0] +16 * q[i+0][j-2][k+1] -1 * q[i+0][j-2][k+2] +8 * q[i+0][j-1][k-2] -128 * q[i+0][j-1][k-1] +240 * q[i+0][j-1][k+0] -128 * q[i+0][j-1][k+1] +8 * q[i+0][j-1][k+2] -8 * q[i+0][j+1][k-2] +128 * q[i+0][j+1][k-1] -240 * q[i+0][j+1][k+0] +128 * q[i+0][j+1][k+1] -8 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -16 * q[i+0][j+2][k-1] +30 * q[i+0][j+2][k+0] -16 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h);}
double d_zzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j+0][k-2] +2 * q[i+0][j+0][k-1] -2 * q[i+0][j+0][k+1] +1 * q[i+0][j+0][k+2])/(2*h*h*h);}
double d_xxxx(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j+0][k+0] -4 * q[i-1][j+0][k+0] +6 * q[i+0][j+0][k+0] -4 * q[i+1][j+0][k+0] +1 * q[i+2][j+0][k+0])/(1*h*h*h*h);}
double d_xxxy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +8 * q[i-2][j-1][k+0] -8 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +2 * q[i-1][j-2][k+0] -16 * q[i-1][j-1][k+0] +16 * q[i-1][j+1][k+0] -2 * q[i-1][j+2][k+0] -2 * q[i+1][j-2][k+0] +16 * q[i+1][j-1][k+0] -16 * q[i+1][j+1][k+0] +2 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -8 * q[i+2][j-1][k+0] +8 * q[i+2][j+1][k+0] -1 * q[i+2][j+2][k+0])/(24*h*h*h*h);}
double d_xxyy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j-2][k+0] -16 * q[i-2][j-1][k+0] +30 * q[i-2][j+0][k+0] -16 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] -16 * q[i-1][j-2][k+0] +256 * q[i-1][j-1][k+0] -480 * q[i-1][j+0][k+0] +256 * q[i-1][j+1][k+0] -16 * q[i-1][j+2][k+0] +30 * q[i+0][j-2][k+0] -480 * q[i+0][j-1][k+0] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+1][k+0] +30 * q[i+0][j+2][k+0] -16 * q[i+1][j-2][k+0] +256 * q[i+1][j-1][k+0] -480 * q[i+1][j+0][k+0] +256 * q[i+1][j+1][k+0] -16 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -16 * q[i+2][j-1][k+0] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h*h);}
double d_xyyy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +2 * q[i-2][j-1][k+0] -2 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +8 * q[i-1][j-2][k+0] -16 * q[i-1][j-1][k+0] +16 * q[i-1][j+1][k+0] -8 * q[i-1][j+2][k+0] -8 * q[i+1][j-2][k+0] +16 * q[i+1][j-1][k+0] -16 * q[i+1][j+1][k+0] +8 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -2 * q[i+2][j-1][k+0] +2 * q[i+2][j+1][k+0] -1 * q[i+2][j+2][k+0])/(24*h*h*h*h);}
double d_yyyy(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j-2][k+0] -4 * q[i+0][j-1][k+0] +6 * q[i+0][j+0][k+0] -4 * q[i+0][j+1][k+0] +1 * q[i+0][j+2][k+0])/(1*h*h*h*h);}
double d_xxxz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +8 * q[i-2][j+0][k-1] -8 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +2 * q[i-1][j+0][k-2] -16 * q[i-1][j+0][k-1] +16 * q[i-1][j+0][k+1] -2 * q[i-1][j+0][k+2] -2 * q[i+1][j+0][k-2] +16 * q[i+1][j+0][k-1] -16 * q[i+1][j+0][k+1] +2 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -8 * q[i+2][j+0][k-1] +8 * q[i+2][j+0][k+1] -1 * q[i+2][j+0][k+2])/(24*h*h*h*h);}
double d_xxyz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +8 * q[i-2][j-2][k-1] -8 * q[i-2][j-2][k+1] +1 * q[i-2][j-2][k+2] +8 * q[i-2][j-1][k-2] -64 * q[i-2][j-1][k-1] +64 * q[i-2][j-1][k+1] -8 * q[i-2][j-1][k+2] -8 * q[i-2][j+1][k-2] +64 * q[i-2][j+1][k-1] -64 * q[i-2][j+1][k+1] +8 * q[i-2][j+1][k+2] +1 * q[i-2][j+2][k-2] -8 * q[i-2][j+2][k-1] +8 * q[i-2][j+2][k+1] -1 * q[i-2][j+2][k+2] +16 * q[i-1][j-2][k-2] -128 * q[i-1][j-2][k-1] +128 * q[i-1][j-2][k+1] -16 * q[i-1][j-2][k+2] -128 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1024 * q[i-1][j-1][k+1] +128 * q[i-1][j-1][k+2] +128 * q[i-1][j+1][k-2] -1024 * q[i-1][j+1][k-1] +1024 * q[i-1][j+1][k+1] -128 * q[i-1][j+1][k+2] -16 * q[i-1][j+2][k-2] +128 * q[i-1][j+2][k-1] -128 * q[i-1][j+2][k+1] +16 * q[i-1][j+2][k+2] -30 * q[i+0][j-2][k-2] +240 * q[i+0][j-2][k-1] -240 * q[i+0][j-2][k+1] +30 * q[i+0][j-2][k+2] +240 * q[i+0][j-1][k-2] -1920 * q[i+0][j-1][k-1] +1920 * q[i+0][j-1][k+1] -240 * q[i+0][j-1][k+2] -240 * q[i+0][j+1][k-2] +1920 * q[i+0][j+1][k-1] -1920 * q[i+0][j+1][k+1] +240 * q[i+0][j+1][k+2] +30 * q[i+0][j+2][k-2] -240 * q[i+0][j+2][k-1] +240 * q[i+0][j+2][k+1] -30 * q[i+0][j+2][k+2] +16 * q[i+1][j-2][k-2] -128 * q[i+1][j-2][k-1] +128 * q[i+1][j-2][k+1] -16 * q[i+1][j-2][k+2] -128 * q[i+1][j-1][k-2] +1024 * q[i+1][j-1][k-1] -1024 * q[i+1][j-1][k+1] +128 * q[i+1][j-1][k+2] +128 * q[i+1][j+1][k-2] -1024 * q[i+1][j+1][k-1] +1024 * q[i+1][j+1][k+1] -128 * q[i+1][j+1][k+2] -16 * q[i+1][j+2][k-2] +128 * q[i+1][j+2][k-1] -128 * q[i+1][j+2][k+1] +16 * q[i+1][j+2][k+2] -1 * q[i+2][j-2][k-2] +8 * q[i+2][j-2][k-1] -8 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] +8 * q[i+2][j-1][k-2] -64 * q[i+2][j-1][k-1] +64 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] -8 * q[i+2][j+1][k-2] +64 * q[i+2][j+1][k-1] -64 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_xyyz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +8 * q[i-2][j-2][k-1] -8 * q[i-2][j-2][k+1] +1 * q[i-2][j-2][k+2] +16 * q[i-2][j-1][k-2] -128 * q[i-2][j-1][k-1] +128 * q[i-2][j-1][k+1] -16 * q[i-2][j-1][k+2] -30 * q[i-2][j+0][k-2] +240 * q[i-2][j+0][k-1] -240 * q[i-2][j+0][k+1] +30 * q[i-2][j+0][k+2] +16 * q[i-2][j+1][k-2] -128 * q[i-2][j+1][k-1] +128 * q[i-2][j+1][k+1] -16 * q[i-2][j+1][k+2] -1 * q[i-2][j+2][k-2] +8 * q[i-2][j+2][k-1] -8 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] +8 * q[i-1][j-2][k-2] -64 * q[i-1][j-2][k-1] +64 * q[i-1][j-2][k+1] -8 * q[i-1][j-2][k+2] -128 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1024 * q[i-1][j-1][k+1] +128 * q[i-1][j-1][k+2] +240 * q[i-1][j+0][k-2] -1920 * q[i-1][j+0][k-1] +1920 * q[i-1][j+0][k+1] -240 * q[i-1][j+0][k+2] -128 * q[i-1][j+1][k-2] +1024 * q[i-1][j+1][k-1] -1024 * q[i-1][j+1][k+1] +128 * q[i-1][j+1][k+2] +8 * q[i-1][j+2][k-2] -64 * q[i-1][j+2][k-1] +64 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] -8 * q[i+1][j-2][k-2] +64 * q[i+1][j-2][k-1] -64 * q[i+1][j-2][k+1] +8 * q[i+1][j-2][k+2] +128 * q[i+1][j-1][k-2] -1024 * q[i+1][j-1][k-1] +1024 * q[i+1][j-1][k+1] -128 * q[i+1][j-1][k+2] -240 * q[i+1][j+0][k-2] +1920 * q[i+1][j+0][k-1] -1920 * q[i+1][j+0][k+1] +240 * q[i+1][j+0][k+2] +128 * q[i+1][j+1][k-2] -1024 * q[i+1][j+1][k-1] +1024 * q[i+1][j+1][k+1] -128 * q[i+1][j+1][k+2] -8 * q[i+1][j+2][k-2] +64 * q[i+1][j+2][k-1] -64 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] +1 * q[i+2][j-2][k-2] -8 * q[i+2][j-2][k-1] +8 * q[i+2][j-2][k+1] -1 * q[i+2][j-2][k+2] -16 * q[i+2][j-1][k-2] +128 * q[i+2][j-1][k-1] -128 * q[i+2][j-1][k+1] +16 * q[i+2][j-1][k+2] +30 * q[i+2][j+0][k-2] -240 * q[i+2][j+0][k-1] +240 * q[i+2][j+0][k+1] -30 * q[i+2][j+0][k+2] -16 * q[i+2][j+1][k-2] +128 * q[i+2][j+1][k-1] -128 * q[i+2][j+1][k+1] +16 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_yyyz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +8 * q[i+0][j-2][k-1] -8 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +2 * q[i+0][j-1][k-2] -16 * q[i+0][j-1][k-1] +16 * q[i+0][j-1][k+1] -2 * q[i+0][j-1][k+2] -2 * q[i+0][j+1][k-2] +16 * q[i+0][j+1][k-1] -16 * q[i+0][j+1][k+1] +2 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -8 * q[i+0][j+2][k-1] +8 * q[i+0][j+2][k+1] -1 * q[i+0][j+2][k+2])/(24*h*h*h*h);}
double d_xxzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i-2][j+0][k-2] -16 * q[i-2][j+0][k-1] +30 * q[i-2][j+0][k+0] -16 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] -16 * q[i-1][j+0][k-2] +256 * q[i-1][j+0][k-1] -480 * q[i-1][j+0][k+0] +256 * q[i-1][j+0][k+1] -16 * q[i-1][j+0][k+2] +30 * q[i+0][j+0][k-2] -480 * q[i+0][j+0][k-1] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] -16 * q[i+1][j+0][k-2] +256 * q[i+1][j+0][k-1] -480 * q[i+1][j+0][k+0] +256 * q[i+1][j+0][k+1] -16 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -16 * q[i+2][j+0][k-1] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h*h);}
double d_xyzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +16 * q[i-2][j-2][k-1] -30 * q[i-2][j-2][k+0] +16 * q[i-2][j-2][k+1] -1 * q[i-2][j-2][k+2] +8 * q[i-2][j-1][k-2] -128 * q[i-2][j-1][k-1] +240 * q[i-2][j-1][k+0] -128 * q[i-2][j-1][k+1] +8 * q[i-2][j-1][k+2] -8 * q[i-2][j+1][k-2] +128 * q[i-2][j+1][k-1] -240 * q[i-2][j+1][k+0] +128 * q[i-2][j+1][k+1] -8 * q[i-2][j+1][k+2] +1 * q[i-2][j+2][k-2] -16 * q[i-2][j+2][k-1] +30 * q[i-2][j+2][k+0] -16 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] +8 * q[i-1][j-2][k-2] -128 * q[i-1][j-2][k-1] +240 * q[i-1][j-2][k+0] -128 * q[i-1][j-2][k+1] +8 * q[i-1][j-2][k+2] -64 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1920 * q[i-1][j-1][k+0] +1024 * q[i-1][j-1][k+1] -64 * q[i-1][j-1][k+2] +64 * q[i-1][j+1][k-2] -1024 * q[i-1][j+1][k-1] +1920 * q[i-1][j+1][k+0] -1024 * q[i-1][j+1][k+1] +64 * q[i-1][j+1][k+2] -8 * q[i-1][j+2][k-2] +128 * q[i-1][j+2][k-1] -240 * q[i-1][j+2][k+0] +128 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] -8 * q[i+1][j-2][k-2] +128 * q[i+1][j-2][k-1] -240 * q[i+1][j-2][k+0] +128 * q[i+1][j-2][k+1] -8 * q[i+1][j-2][k+2] +64 * q[i+1][j-1][k-2] -1024 * q[i+1][j-1][k-1] +1920 * q[i+1][j-1][k+0] -1024 * q[i+1][j-1][k+1] +64 * q[i+1][j-1][k+2] -64 * q[i+1][j+1][k-2] +1024 * q[i+1][j+1][k-1] -1920 * q[i+1][j+1][k+0] +1024 * q[i+1][j+1][k+1] -64 * q[i+1][j+1][k+2] +8 * q[i+1][j+2][k-2] -128 * q[i+1][j+2][k-1] +240 * q[i+1][j+2][k+0] -128 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] +1 * q[i+2][j-2][k-2] -16 * q[i+2][j-2][k-1] +30 * q[i+2][j-2][k+0] -16 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] -8 * q[i+2][j-1][k-2] +128 * q[i+2][j-1][k-1] -240 * q[i+2][j-1][k+0] +128 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] +8 * q[i+2][j+1][k-2] -128 * q[i+2][j+1][k-1] +240 * q[i+2][j+1][k+0] -128 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] -1 * q[i+2][j+2][k-2] +16 * q[i+2][j+2][k-1] -30 * q[i+2][j+2][k+0] +16 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_yyzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j-2][k-2] -16 * q[i+0][j-2][k-1] +30 * q[i+0][j-2][k+0] -16 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] -16 * q[i+0][j-1][k-2] +256 * q[i+0][j-1][k-1] -480 * q[i+0][j-1][k+0] +256 * q[i+0][j-1][k+1] -16 * q[i+0][j-1][k+2] +30 * q[i+0][j+0][k-2] -480 * q[i+0][j+0][k-1] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] -16 * q[i+0][j+1][k-2] +256 * q[i+0][j+1][k-1] -480 * q[i+0][j+1][k+0] +256 * q[i+0][j+1][k+1] -16 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -16 * q[i+0][j+2][k-1] +30 * q[i+0][j+2][k+0] -16 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h*h);}
double d_xzzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +2 * q[i-2][j+0][k-1] -2 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +8 * q[i-1][j+0][k-2] -16 * q[i-1][j+0][k-1] +16 * q[i-1][j+0][k+1] -8 * q[i-1][j+0][k+2] -8 * q[i+1][j+0][k-2] +16 * q[i+1][j+0][k-1] -16 * q[i+1][j+0][k+1] +8 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -2 * q[i+2][j+0][k-1] +2 * q[i+2][j+0][k+1] -1 * q[i+2][j+0][k+2])/(24*h*h*h*h);}
double d_yzzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +2 * q[i+0][j-2][k-1] -2 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +8 * q[i+0][j-1][k-2] -16 * q[i+0][j-1][k-1] +16 * q[i+0][j-1][k+1] -8 * q[i+0][j-1][k+2] -8 * q[i+0][j+1][k-2] +16 * q[i+0][j+1][k-1] -16 * q[i+0][j+1][k+1] +8 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -2 * q[i+0][j+2][k-1] +2 * q[i+0][j+2][k+1] -1 * q[i+0][j+2][k+2])/(24*h*h*h*h);}
double d_zzzz(int i,int j,int k,double q[MX*NX+2*Ns][MY*NY+2*Ns][MZ*NZ+2*Ns]) { return (+1 * q[i+0][j+0][k-2] -4 * q[i+0][j+0][k-1] +6 * q[i+0][j+0][k+0] -4 * q[i+0][j+0][k+1] +1 * q[i+0][j+0][k+2])/(1*h*h*h*h);}

void init(navi &n, floors<double> &s) {
  unique_ptr<tmp_buf> tmp(new tmp_buf);

  // naviの初期化
  n.time_step = 0;
  n.lower_x = 2*NX;
  n.lower_y = 2*NY;
  n.lower_z = 2*NZ;
  n.upper_x = (MX+2)*NX;
  n.upper_y = (MX+2)*NY;
  n.upper_z = (MX+2)*NZ;
  n.offset_x = -2*NX;
  n.offset_y = -2*NY;
  n.offset_z = -2*NZ;

#if defined(VORTEX)
  // 初期化 (一様等方乱流)
  setup_vortex();
#else
  // 初期化 (一様等方乱流)
  setup_turbulence(*tmp);
#endif

  // predictor, half-correctorの構築
  for(int ix = Ns; ix < LX+Ns; ix++) {
    int x = ix - Ns + 2*NX;
    for(int iy = Ns; iy < LY+Ns; iy++) {
      int y = iy - Ns + 2*NY;
      for(int iz = Ns; iz < LZ+Ns; iz++) {
        int z = iz - Ns + 2*NZ;
        double r = tmp->rc[ix][iy][iz];
        double u = tmp->uc[ix][iy][iz];
        double v = tmp->vc[ix][iy][iz];
        double w = tmp->wc[ix][iy][iz];
        double p = tmp->pc[ix][iy][iz];

        double r_x = d_x(ix,iy,iz,tmp->rc);
        double u_x = d_x(ix,iy,iz,tmp->uc);
        double v_x = d_x(ix,iy,iz,tmp->vc);
        double w_x = d_x(ix,iy,iz,tmp->wc);
        double p_x = d_x(ix,iy,iz,tmp->pc);
        double r_y = d_y(ix,iy,iz,tmp->rc);
        double u_y = d_y(ix,iy,iz,tmp->uc);
        double v_y = d_y(ix,iy,iz,tmp->vc);
        double w_y = d_y(ix,iy,iz,tmp->wc);
        double p_y = d_y(ix,iy,iz,tmp->pc);
        double r_z = d_z(ix,iy,iz,tmp->rc);
        double u_z = d_z(ix,iy,iz,tmp->uc);
        double v_z = d_z(ix,iy,iz,tmp->vc);
        double w_z = d_z(ix,iy,iz,tmp->wc);
        double p_z = d_z(ix,iy,iz,tmp->pc);

        double r_xx = d_xx(ix,iy,iz,tmp->rc);
        double u_xx = d_xx(ix,iy,iz,tmp->uc);
        double v_xx = d_xx(ix,iy,iz,tmp->vc);
        double w_xx = d_xx(ix,iy,iz,tmp->wc);
        double p_xx = d_xx(ix,iy,iz,tmp->pc);
        double r_yy = d_yy(ix,iy,iz,tmp->rc);
        double u_yy = d_yy(ix,iy,iz,tmp->uc);
        double v_yy = d_yy(ix,iy,iz,tmp->vc);
        double w_yy = d_yy(ix,iy,iz,tmp->wc);
        double p_yy = d_yy(ix,iy,iz,tmp->pc);
        double r_zz = d_zz(ix,iy,iz,tmp->rc);
        double u_zz = d_zz(ix,iy,iz,tmp->uc);
        double v_zz = d_zz(ix,iy,iz,tmp->vc);
        double w_zz = d_zz(ix,iy,iz,tmp->wc);
        double p_zz = d_zz(ix,iy,iz,tmp->pc);
        double r_xy = d_xy(ix,iy,iz,tmp->rc);
        double u_xy = d_xy(ix,iy,iz,tmp->uc);
        double v_xy = d_xy(ix,iy,iz,tmp->vc);
        double w_xy = d_xy(ix,iy,iz,tmp->wc);
        double p_xy = d_xy(ix,iy,iz,tmp->pc);
        double r_yz = d_yz(ix,iy,iz,tmp->rc);
        double u_yz = d_yz(ix,iy,iz,tmp->uc);
        double v_yz = d_yz(ix,iy,iz,tmp->vc);
        double w_yz = d_yz(ix,iy,iz,tmp->wc);
        double p_yz = d_yz(ix,iy,iz,tmp->pc);
        double r_xz = d_xz(ix,iy,iz,tmp->rc);
        double u_xz = d_xz(ix,iy,iz,tmp->uc);
        double v_xz = d_xz(ix,iy,iz,tmp->vc);
        double w_xz = d_xz(ix,iy,iz,tmp->wc);
        double p_xz = d_xz(ix,iy,iz,tmp->pc);

        double u_xxx = d_xxx(ix,iy,iz,tmp->uc);
        double v_xxx = d_xxx(ix,iy,iz,tmp->vc);
        double w_xxx = d_xxx(ix,iy,iz,tmp->wc);
        double p_xxx = d_xxx(ix,iy,iz,tmp->pc);
        double u_yyy = d_yyy(ix,iy,iz,tmp->uc);
        double v_yyy = d_yyy(ix,iy,iz,tmp->vc);
        double w_yyy = d_yyy(ix,iy,iz,tmp->wc);
        double p_yyy = d_yyy(ix,iy,iz,tmp->pc);
        double u_zzz = d_zzz(ix,iy,iz,tmp->uc);
        double v_zzz = d_zzz(ix,iy,iz,tmp->vc);
        double w_zzz = d_zzz(ix,iy,iz,tmp->wc);
        double p_zzz = d_zzz(ix,iy,iz,tmp->pc);
        double u_xxy = d_xxy(ix,iy,iz,tmp->uc);
        double v_xxy = d_xxy(ix,iy,iz,tmp->vc);
        double w_xxy = d_xxy(ix,iy,iz,tmp->wc);
        double p_xxy = d_xxy(ix,iy,iz,tmp->pc);
        double u_xxz = d_xxz(ix,iy,iz,tmp->uc);
        double v_xxz = d_xxz(ix,iy,iz,tmp->vc);
        double w_xxz = d_xxz(ix,iy,iz,tmp->wc);
        double p_xxz = d_xxz(ix,iy,iz,tmp->pc);
        double u_xyy = d_xyy(ix,iy,iz,tmp->uc);
        double v_xyy = d_xyy(ix,iy,iz,tmp->vc);
        double w_xyy = d_xyy(ix,iy,iz,tmp->wc);
        double p_xyy = d_xyy(ix,iy,iz,tmp->pc);
        double u_xzz = d_xzz(ix,iy,iz,tmp->uc);
        double v_xzz = d_xzz(ix,iy,iz,tmp->vc);
        double w_xzz = d_xzz(ix,iy,iz,tmp->wc);
        double p_xzz = d_xzz(ix,iy,iz,tmp->pc);
        double u_yyz = d_yyz(ix,iy,iz,tmp->uc);
        double v_yyz = d_yyz(ix,iy,iz,tmp->vc);
        double w_yyz = d_yyz(ix,iy,iz,tmp->wc);
        double p_yyz = d_yyz(ix,iy,iz,tmp->pc);
        double u_yzz = d_yzz(ix,iy,iz,tmp->uc);
        double v_yzz = d_yzz(ix,iy,iz,tmp->vc);
        double w_yzz = d_yzz(ix,iy,iz,tmp->wc);
        double p_yzz = d_yzz(ix,iy,iz,tmp->pc);
        double u_xyz = d_xyz(ix,iy,iz,tmp->uc);
        double v_xyz = d_xyz(ix,iy,iz,tmp->vc);
        double w_xyz = d_xyz(ix,iy,iz,tmp->wc);

        double u_xxxx = d_xxxx(ix,iy,iz,tmp->uc);
        double u_xxxy = d_xxxy(ix,iy,iz,tmp->uc);
        double u_xxyy = d_xxyy(ix,iy,iz,tmp->uc);
        double u_xyyy = d_xyyy(ix,iy,iz,tmp->uc);
        double u_yyyy = d_yyyy(ix,iy,iz,tmp->uc);
        double u_xxxz = d_xxxz(ix,iy,iz,tmp->uc);
        double u_xyyz = d_xyyz(ix,iy,iz,tmp->uc);
        double u_xxzz = d_xxzz(ix,iy,iz,tmp->uc);
        double u_xyzz = d_xyzz(ix,iy,iz,tmp->uc);
        double u_yyzz = d_yyzz(ix,iy,iz,tmp->uc);
        double u_xzzz = d_xzzz(ix,iy,iz,tmp->uc);
        double u_zzzz = d_zzzz(ix,iy,iz,tmp->uc);

        double v_xxxx = d_xxxx(ix,iy,iz,tmp->vc);
        double v_xxxy = d_xxxy(ix,iy,iz,tmp->vc);
        double v_xxyy = d_xxyy(ix,iy,iz,tmp->vc);
        double v_xyyy = d_xyyy(ix,iy,iz,tmp->vc);
        double v_yyyy = d_yyyy(ix,iy,iz,tmp->vc);
        double v_xxyz = d_xxyz(ix,iy,iz,tmp->vc);
        double v_yyyz = d_yyyz(ix,iy,iz,tmp->vc);
        double v_xxzz = d_xxzz(ix,iy,iz,tmp->vc);
        double v_xyzz = d_xyzz(ix,iy,iz,tmp->vc);
        double v_yyzz = d_yyzz(ix,iy,iz,tmp->vc);
        double v_yzzz = d_yzzz(ix,iy,iz,tmp->vc);
        double v_zzzz = d_zzzz(ix,iy,iz,tmp->vc);

        double w_xxxx = d_xxxx(ix,iy,iz,tmp->wc);
        double w_xxyy = d_xxyy(ix,iy,iz,tmp->wc);
        double w_yyyy = d_yyyy(ix,iy,iz,tmp->wc);
        double w_xxxz = d_xxxz(ix,iy,iz,tmp->wc);
        double w_xxyz = d_xxyz(ix,iy,iz,tmp->wc);
        double w_xyyz = d_xyyz(ix,iy,iz,tmp->wc);
        double w_yyyz = d_yyyz(ix,iy,iz,tmp->wc);
        double w_xxzz = d_xxzz(ix,iy,iz,tmp->wc);
        double w_yyzz = d_yyzz(ix,iy,iz,tmp->wc);
        double w_xzzz = d_xzzz(ix,iy,iz,tmp->wc);
        double w_yzzz = d_yzzz(ix,iy,iz,tmp->wc);
        double w_zzzz = d_zzzz(ix,iy,iz,tmp->wc);

        double vis1 = 4.0/3.0*u_xx + 1.0/3.0*v_xy + 1.0/3.0*w_xz + u_yy + u_zz;
        double vis2 = 1.0/3.0*u_xy + 4.0/3.0*v_yy + 1.0/3.0*w_yz + v_xx + v_zz;
        double vis3 = 1.0/3.0*u_xz + 1.0/3.0*v_yz + 4.0/3.0*w_zz + w_xx + w_yy;
        double vis1_x = 4.0/3.0*u_xxx + 1.0/3.0*v_xxy + 1.0/3.0*w_xxz + u_xyy + u_xzz;
        double vis2_x = 1.0/3.0*u_xxy + 4.0/3.0*v_xyy + 1.0/3.0*w_xyz + v_xzz + v_xxx;
        double vis3_x = 1.0/3.0*u_xxz + 1.0/3.0*v_xyz + 4.0/3.0*w_xzz + w_xyy + w_xxx;
        double vis1_y = 4.0/3.0*u_xxy + 1.0/3.0*v_xyy + 1.0/3.0*w_xyz + u_yzz + u_yyy;
        double vis2_y = 1.0/3.0*u_xyy + 4.0/3.0*v_yyy + 1.0/3.0*w_yyz + v_xxy + v_yzz;
        double vis3_y = 1.0/3.0*u_xyz + 1.0/3.0*v_yyz + 4.0/3.0*w_yzz + w_xxy + w_yyy;
        double vis1_z = 4.0/3.0*u_xxz + 1.0/3.0*v_xyz + 1.0/3.0*w_xzz + u_yyz + u_zzz;
        double vis2_z = 1.0/3.0*u_xyz + 4.0/3.0*v_yyz + 1.0/3.0*w_yzz + v_xxz + v_zzz;
        double vis3_z = 1.0/3.0*u_xzz + 1.0/3.0*v_yzz + 4.0/3.0*w_zzz + w_xxz + w_yyz;
        double vis1_xx = 4.0/3.0*u_xxxx + 1.0/3.0*v_xxxy + 1.0/3.0*w_xxxz + u_xxyy + u_xxzz;
        double vis2_xx = 1.0/3.0*u_xxxy + 4.0/3.0*v_xxyy + 1.0/3.0*w_xxyz + v_xxzz + v_xxxx;
        double vis3_xx = 1.0/3.0*u_xxxz + 1.0/3.0*v_xxyz + 4.0/3.0*w_xxzz + w_xxyy + w_xxxx;
        double vis1_yy = 4.0/3.0*u_xxyy + 1.0/3.0*v_xyyy + 1.0/3.0*w_xyyz + u_yyzz + u_yyyy;
        double vis2_yy = 1.0/3.0*u_xyyy + 4.0/3.0*v_yyyy + 1.0/3.0*w_yyyz + v_xxyy + v_yyzz;
        double vis3_yy = 1.0/3.0*u_xyyz + 1.0/3.0*v_yyyz + 4.0/3.0*w_yyzz + w_xxyy + w_yyyy;
        double vis1_zz = 4.0/3.0*u_xxzz + 1.0/3.0*v_xyzz + 1.0/3.0*w_xzzz + u_yyzz + u_zzzz;
        double vis2_zz = 1.0/3.0*u_xyzz + 4.0/3.0*v_yyzz + 1.0/3.0*w_yzzz + v_xxzz + v_zzzz;
        double vis3_zz = 1.0/3.0*u_xzzz + 1.0/3.0*v_yzzz + 4.0/3.0*w_zzzz + w_xxzz + w_yyzz;
        double vis1_xy = 4.0/3.0*u_xxxy + 1.0/3.0*v_xxyy + 1.0/3.0*w_xxyz + u_xyzz + u_xyyy;
        double vis2_xy = 1.0/3.0*u_xxyy + 4.0/3.0*v_xyyy + 1.0/3.0*w_xyyz + v_xyzz + v_xxxy;
        double vis2_yz = 1.0/3.0*u_xyyz + 4.0/3.0*v_yyyz + 1.0/3.0*w_yyzz + v_xxyz + v_yzzz;
        double vis3_yz = 1.0/3.0*u_xyzz + 1.0/3.0*v_yyzz + 4.0/3.0*w_yzzz + w_xxyz + w_yyyz;
        double vis1_xz = 4.0/3.0*u_xxxz + 1.0/3.0*v_xxyz + 1.0/3.0*w_xxzz + u_xyyz + u_xzzz;
        double vis3_xz = 1.0/3.0*u_xxzz + 1.0/3.0*v_xyzz + 4.0/3.0*w_xzzz + w_xyyz + w_xxxz;

        double r_t = -r*u_x - r*v_y - r*w_z - r_x*u - r_y*v - r_z*w;
        double u_t = c*pow(r,-1.0)*vis1 - p_x*pow(r,-1.0) - u*u_x - u_y*v - u_z*w;
        double v_t = c*pow(r,-1.0)*vis2 - p_y*pow(r,-1.0) - u*v_x - v*v_y - v_z*w;
        double w_t = c*pow(r,-1.0)*vis3 - p_z*pow(r,-1.0) - u*w_x - v*w_y - w*w_z;
        double p_t = -c2*u*vis1 - c2*v*vis2 - c2*vis3*w - gm*p*u_x - gm*p*v_y - gm*p*w_z - p_x*u - p_y*v - p_z*w;

        double r_tx = -r*u_xx - r*v_xy - r*w_xz - 2.0*r_x*u_x - r_x*v_y - r_x*w_z - r_xy*v - r_xz*w - r_xx*u - r_y*v_x - r_z*w_x;
        double u_tx = -c*pow(r,-2.0)*r_x*vis1 + c*pow(r,-1.0)*vis1_x + p_x*pow(r,-2.0)*r_x - p_xx*pow(r,-1.0) - u*u_xx - pow(u_x,2.0) - u_xy*v - u_xz*w - u_y*v_x - u_z*w_x;
        double v_tx = -c*pow(r,-2.0)*r_x*vis2 + c*pow(r,-1.0)*vis2_x - p_xy*pow(r,-1.0) + p_y*pow(r,-2.0)*r_x - u*v_xx - u_x*v_x - v*v_xy - v_x*v_y - v_xz*w - v_z*w_x;
        double w_tx = -c*pow(r,-2.0)*r_x*vis3 + c*pow(r,-1.0)*vis3_x - p_xz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_x - u*w_xx - u_x*w_x - v*w_xy - v_x*w_y - w*w_xz - w_x*w_z;
        double p_tx = -c2*u*vis1_x - c2*u_x*vis1 - c2*v*vis2_x - c2*v_x*vis2 - c2*vis3*w_x - c2*vis3_x*w - gm*p*u_xx - gm*p*v_xy - gm*p*w_xz - gm*p_x*u_x - gm*p_x*v_y - gm*p_x*w_z - p_x*u_x - p_xy*v - p_xz*w - p_xx*u - p_y*v_x - p_z*w_x;

        double r_ty = -r*u_xy - r*v_yy - r*w_yz - r_x*u_y - r_xy*u - r_y*u_x - 2.0*r_y*v_y - r_y*w_z - r_yz*w - r_yy*v - r_z*w_y;
        double u_ty = -c*pow(r,-2.0)*r_y*vis1 + c*pow(r,-1.0)*vis1_y + p_x*pow(r,-2.0)*r_y - p_xy*pow(r,-1.0) - u*u_xy - u_x*u_y - u_y*v_y - u_yz*w - u_yy*v - u_z*w_y;
        double v_ty = -c*pow(r,-2.0)*r_y*vis2 + c*pow(r,-1.0)*vis2_y + p_y*pow(r,-2.0)*r_y - p_yy*pow(r,-1.0) - u*v_xy - u_y*v_x - v*v_yy - pow(v_y,2.0) - v_yz*w - v_z*w_y;
        double w_ty = -c*pow(r,-2.0)*r_y*vis3 + c*pow(r,-1.0)*vis3_y - p_yz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_y - u*w_xy - u_y*w_x - v*w_yy - v_y*w_y - w*w_yz - w_y*w_z;
        double p_ty = -c2*u*vis1_y - c2*u_y*vis1 - c2*v*vis2_y - c2*v_y*vis2 - c2*vis3*w_y - c2*vis3_y*w - gm*p*u_xy - gm*p*v_yy - gm*p*w_yz - gm*p_y*u_x - gm*p_y*v_y - gm*p_y*w_z - p_x*u_y - p_xy*u - p_y*v_y - p_yz*w - p_yy*v - p_z*w_y;

        double r_tz = -r*u_xz - r*v_yz - r*w_zz - r_x*u_z - r_xz*u - r_y*v_z - r_yz*v - r_z*u_x - r_z*v_y - 2.0*r_z*w_z - r_zz*w;
        double u_tz = -c*pow(r,-2.0)*r_z*vis1 + c*pow(r,-1.0)*vis1_z + p_x*pow(r,-2.0)*r_z - p_xz*pow(r,-1.0) - u*u_xz - u_x*u_z - u_y*v_z - u_yz*v - u_z*w_z - u_zz*w;
        double v_tz = -c*pow(r,-2.0)*r_z*vis2 + c*pow(r,-1.0)*vis2_z + p_y*pow(r,-2.0)*r_z - p_yz*pow(r,-1.0) - u*v_xz - u_z*v_x - v*v_yz - v_y*v_z - v_z*w_z - v_zz*w;
        double w_tz = -c*pow(r,-2.0)*r_z*vis3 + c*pow(r,-1.0)*vis3_z + p_z*pow(r,-2.0)*r_z - p_zz*pow(r,-1.0) - u*w_xz - u_z*w_x - v*w_yz - v_z*w_y - w*w_zz - pow(w_z,2.0);
        double p_tz = -c2*u*vis1_z - c2*u_z*vis1 - c2*v*vis2_z - c2*v_z*vis2 - c2*vis3*w_z - c2*vis3_z*w - gm*p*u_xz - gm*p*v_yz - gm*p*w_zz - gm*p_z*u_x - gm*p_z*v_y - gm*p_z*w_z - p_x*u_z - p_xz*u - p_y*v_z - p_yz*v - p_z*w_z - p_zz*w;

        double u_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_x*vis1_x - c*pow(r,-2.0)*r_xx*vis1 + c*pow(r,-1.0)*vis1_xx - 2.0*p_x*pow(r,-3.0)*pow(r_x,2.0) + p_x*pow(r,-2.0)*r_xx + 2.0*p_xx*pow(r,-2.0)*r_x - p_xxx*pow(r,-1.0) - u*u_xxx - 3.0*u_x*u_xx - 2.0*u_xy*v_x - 2.0*u_xz*w_x - u_xxy*v - u_xxz*w - u_y*v_xx - u_z*w_xx;
        double v_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_x*vis2_x - c*pow(r,-2.0)*r_xx*vis2 + c*pow(r,-1.0)*vis2_xx + 2.0*p_xy*pow(r,-2.0)*r_x - p_xxy*pow(r,-1.0) - 2.0*p_y*pow(r,-3.0)*pow(r_x,2.0) + p_y*pow(r,-2.0)*r_xx - u*v_xxx - 2.0*u_x*v_xx - u_xx*v_x - v*v_xxy - 2.0*v_x*v_xy - 2.0*v_xz*w_x - v_xx*v_y - v_xxz*w - v_z*w_xx;
        double w_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_x*vis3_x - c*pow(r,-2.0)*r_xx*vis3 + c*pow(r,-1.0)*vis3_xx + 2.0*p_xz*pow(r,-2.0)*r_x - p_xxz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*pow(r_x,2.0) + p_z*pow(r,-2.0)*r_xx - u*w_xxx - 2.0*u_x*w_xx - u_xx*w_x - v*w_xxy - 2.0*v_x*w_xy - v_xx*w_y - w*w_xxz - 2.0*w_x*w_xz - w_xx*w_z;

        double u_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_y*vis1_y - c*pow(r,-2.0)*r_yy*vis1 + c*pow(r,-1.0)*vis1_yy - 2.0*p_x*pow(r,-3.0)*pow(r_y,2.0) + p_x*pow(r,-2.0)*r_yy + 2.0*p_xy*pow(r,-2.0)*r_y - p_xyy*pow(r,-1.0) - u*u_xyy - u_x*u_yy - 2.0*u_xy*u_y - u_y*v_yy - 2.0*u_yz*w_y - 2.0*u_yy*v_y - u_yyz*w - u_yyy*v - u_z*w_yy;
        double v_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_y*vis2_y - c*pow(r,-2.0)*r_yy*vis2 + c*pow(r,-1.0)*vis2_yy - 2.0*p_y*pow(r,-3.0)*pow(r_y,2.0) + p_y*pow(r,-2.0)*r_yy + 2.0*p_yy*pow(r,-2.0)*r_y - p_yyy*pow(r,-1.0) - u*v_xyy - 2.0*u_y*v_xy - u_yy*v_x - v*v_yyy - 3.0*v_y*v_yy - 2.0*v_yz*w_y - v_yyz*w - v_z*w_yy;
        double w_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_y*vis3_y - c*pow(r,-2.0)*r_yy*vis3 + c*pow(r,-1.0)*vis3_yy + 2.0*p_yz*pow(r,-2.0)*r_y - p_yyz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*pow(r_y,2.0) + p_z*pow(r,-2.0)*r_yy - u*w_xyy - 2.0*u_y*w_xy - u_yy*w_x - v*w_yyy - 2.0*v_y*w_yy - v_yy*w_y - w*w_yyz - 2.0*w_y*w_yz - w_yy*w_z;

        double u_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_z*vis1_z - c*pow(r,-2.0)*r_zz*vis1 + c*pow(r,-1.0)*vis1_zz - 2.0*p_x*pow(r,-3.0)*pow(r_z,2.0) + p_x*pow(r,-2.0)*r_zz + 2.0*p_xz*pow(r,-2.0)*r_z - p_xzz*pow(r,-1.0) - u*u_xzz - u_x*u_zz - 2.0*u_xz*u_z - u_y*v_zz - 2.0*u_yz*v_z - u_yzz*v - u_z*w_zz - 2.0*u_zz*w_z - u_zzz*w;
        double v_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_z*vis2_z - c*pow(r,-2.0)*r_zz*vis2 + c*pow(r,-1.0)*vis2_zz - 2.0*p_y*pow(r,-3.0)*pow(r_z,2.0) + p_y*pow(r,-2.0)*r_zz + 2.0*p_yz*pow(r,-2.0)*r_z - p_yzz*pow(r,-1.0) - u*v_xzz - 2.0*u_z*v_xz - u_zz*v_x - v*v_yzz - v_y*v_zz - 2.0*v_yz*v_z - v_z*w_zz - 2.0*v_zz*w_z - v_zzz*w;
        double w_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_z*vis3_z - c*pow(r,-2.0)*r_zz*vis3 + c*pow(r,-1.0)*vis3_zz - 2.0*p_z*pow(r,-3.0)*pow(r_z,2.0) + p_z*pow(r,-2.0)*r_zz + 2.0*p_zz*pow(r,-2.0)*r_z - p_zzz*pow(r,-1.0) - u*w_xzz - 2.0*u_z*w_xz - u_zz*w_x - v*w_yzz - 2.0*v_z*w_yz - v_zz*w_y - w*w_zzz - 3.0*w_z*w_zz;

        double u_txy = 2.0*c*pow(r,-3.0)*r_x*r_y*vis1 - c*pow(r,-2.0)*r_x*vis1_y - c*pow(r,-2.0)*r_xy*vis1 - c*pow(r,-2.0)*r_y*vis1_x + c*pow(r,-1.0)*vis1_xy - 2.0*p_x*pow(r,-3.0)*r_x*r_y + p_x*pow(r,-2.0)*r_xy + p_xy*pow(r,-2.0)*r_x + p_xx*pow(r,-2.0)*r_y - p_xxy*pow(r,-1.0) - u*u_xxy - 2.0*u_x*u_xy - u_xy*v_y - u_xyz*w - u_xyy*v - u_xz*w_y - u_xx*u_y - u_y*v_xy - u_yz*w_x - u_yy*v_x - u_z*w_xy;
        double v_txy = 2.0*c*pow(r,-3.0)*r_x*r_y*vis2 - c*pow(r,-2.0)*r_x*vis2_y - c*pow(r,-2.0)*r_xy*vis2 - c*pow(r,-2.0)*r_y*vis2_x + c*pow(r,-1.0)*vis2_xy + p_xy*pow(r,-2.0)*r_y - p_xyy*pow(r,-1.0) - 2.0*p_y*pow(r,-3.0)*r_x*r_y + p_y*pow(r,-2.0)*r_xy + p_yy*pow(r,-2.0)*r_x - u*v_xxy - u_x*v_xy - u_xy*v_x - u_y*v_xx - v*v_xyy - v_x*v_yy - 2.0*v_xy*v_y - v_xyz*w - v_xz*w_y - v_yz*w_x - v_z*w_xy;

        double v_tyz = 2.0*c*pow(r,-3.0)*r_y*r_z*vis2 - c*pow(r,-2.0)*r_y*vis2_z - c*pow(r,-2.0)*r_yz*vis2 - c*pow(r,-2.0)*r_z*vis2_y + c*pow(r,-1.0)*vis2_yz - 2.0*p_y*pow(r,-3.0)*r_y*r_z + p_y*pow(r,-2.0)*r_yz + p_yz*pow(r,-2.0)*r_y + p_yy*pow(r,-2.0)*r_z - p_yyz*pow(r,-1.0) - u*v_xyz - u_y*v_xz - u_yz*v_x - u_z*v_xy - v*v_yyz - 2.0*v_y*v_yz - v_yz*w_z - v_yzz*w - v_yy*v_z - v_z*w_yz - v_zz*w_y;
        double w_tyz = 2.0*c*pow(r,-3.0)*r_y*r_z*vis3 - c*pow(r,-2.0)*r_y*vis3_z - c*pow(r,-2.0)*r_yz*vis3 - c*pow(r,-2.0)*r_z*vis3_y + c*pow(r,-1.0)*vis3_yz + p_yz*pow(r,-2.0)*r_z - p_yzz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*r_y*r_z + p_z*pow(r,-2.0)*r_yz + p_zz*pow(r,-2.0)*r_y - u*w_xyz - u_y*w_xz - u_yz*w_x - u_z*w_xy - v*w_yyz - v_y*w_yz - v_yz*w_y - v_z*w_yy - w*w_yzz - w_y*w_zz - 2.0*w_yz*w_z;

        double u_txz = 2.0*c*pow(r,-3.0)*r_x*r_z*vis1 - c*pow(r,-2.0)*r_x*vis1_z - c*pow(r,-2.0)*r_xz*vis1 - c*pow(r,-2.0)*r_z*vis1_x + c*pow(r,-1.0)*vis1_xz - 2.0*p_x*pow(r,-3.0)*r_x*r_z + p_x*pow(r,-2.0)*r_xz + p_xz*pow(r,-2.0)*r_x + p_xx*pow(r,-2.0)*r_z - p_xxz*pow(r,-1.0) - u*u_xxz - 2.0*u_x*u_xz - u_xy*v_z - u_xyz*v - u_xz*w_z - u_xzz*w - u_xx*u_z - u_y*v_xz - u_yz*v_x - u_z*w_xz - u_zz*w_x;
        double w_txz = 2.0*c*pow(r,-3.0)*r_x*r_z*vis3 - c*pow(r,-2.0)*r_x*vis3_z - c*pow(r,-2.0)*r_xz*vis3 - c*pow(r,-2.0)*r_z*vis3_x + c*pow(r,-1.0)*vis3_xz + p_xz*pow(r,-2.0)*r_z - p_xzz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*r_x*r_z + p_z*pow(r,-2.0)*r_xz + p_zz*pow(r,-2.0)*r_x - u*w_xxz - u_x*w_xz - u_xz*w_x - u_z*w_xx - v*w_xyz - v_x*w_yz - v_xz*w_y - v_z*w_xy - w*w_xzz - w_x*w_zz - 2.0*w_xz*w_z;

        double vis1_t = 4.0/3.0*u_txx + 1.0/3.0*v_txy + 1.0/3.0*w_txz + u_tyy + u_tzz;
        double vis2_t = 1.0/3.0*u_txy + 4.0/3.0*v_tyy + 1.0/3.0*w_tyz + v_txx + v_tzz;
        double vis3_t = 1.0/3.0*u_txz + 1.0/3.0*v_tyz + 4.0/3.0*w_tzz + w_txx + w_tyy;

        double r_tt = -r*u_tx - r*v_ty - r*w_tz - r_t*u_x - r_t*v_y - r_t*w_z - r_tx*u - r_ty*v - r_tz*w - r_x*u_t - r_y*v_t - r_z*w_t;
        double u_tt = -c*pow(r,-2.0)*r_t*vis1 + c*pow(r,-1.0)*vis1_t - p_tx*pow(r,-1.0) + p_x*pow(r,-2.0)*r_t - u*u_tx - u_t*u_x - u_ty*v - u_tz*w - u_y*v_t - u_z*w_t;
        double v_tt = -c*pow(r,-2.0)*r_t*vis2 + c*pow(r,-1.0)*vis2_t - p_ty*pow(r,-1.0) + p_y*pow(r,-2.0)*r_t - u*v_tx - u_t*v_x - v*v_ty - v_t*v_y - v_tz*w - v_z*w_t;
        double w_tt = -c*pow(r,-2.0)*r_t*vis3 + c*pow(r,-1.0)*vis3_t - p_tz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_t - u*w_tx - u_t*w_x - v*w_ty - v_t*w_y - w*w_tz - w_t*w_z;
        double p_tt = -c2*u*vis1_t - c2*u_t*vis1 - c2*v*vis2_t - c2*v_t*vis2 - c2*vis3*w_t - c2*vis3_t*w - gm*p*u_tx - gm*p*v_ty - gm*p*w_tz - gm*p_t*u_x - gm*p_t*v_y - gm*p_t*w_z - p_tx*u - p_ty*v - p_tz*w - p_x*u_t - p_y*v_t - p_z*w_t;

        s.p[x][y][z].r = r + dt*r_t + dt*dt*r_tt/2;
        s.p[x][y][z].u = u + dt*u_t + dt*dt*u_tt/2;
        s.p[x][y][z].v = v + dt*v_t + dt*dt*v_tt/2;
        s.p[x][y][z].w = w + dt*w_t + dt*dt*w_tt/2;
        s.p[x][y][z].p = p + dt*p_t + dt*dt*p_tt/2;

        s.h[x][y][z].r = r + dt*r_t/2 - dt*dt*r_tt/12;
        s.h[x][y][z].u = u + dt*u_t/2 - dt*dt*u_tt/12;
        s.h[x][y][z].v = v + dt*v_t/2 - dt*dt*v_tt/12;
        s.h[x][y][z].w = w + dt*w_t/2 - dt*dt*w_tt/12;
        s.h[x][y][z].p = p + dt*p_t/2 - dt*dt*p_tt/12;
      }
    }
  }
}

