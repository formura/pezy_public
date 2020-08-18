#include <memory>
#include "config.h"

using namespace std;

// 6次ローパスフィルター
// 袖3を必要とする

// フィルターの袖サイズ
#define Ms 3

double filter(int ix,int iy,int iz,double q[MX*NX+2*Ms][MY*NY+2*Ms][MZ*NZ+2*Ms]) {
  double a0 = 11.0/16.0;
  double a1 = 15.0/64.0;
  double a2 = -3.0/32.0;
  double a3 = 1.0/64.0;

  return a0*a0*a0*(q[ix+0][iy+0][iz+0])
       + a1*a0*a0*(q[ix-1][iy+0][iz+0] + q[ix+0][iy-1][iz+0] + q[ix+0][iy+0][iz-1] + q[ix+0][iy+0][iz+1] + q[ix+0][iy+1][iz+0] + q[ix+1][iy+0][iz+0])
       + a1*a1*a0*(q[ix-1][iy-1][iz+0] + q[ix-1][iy+0][iz-1] + q[ix-1][iy+0][iz+1] + q[ix-1][iy+1][iz+0] + q[ix+0][iy-1][iz-1] + q[ix+0][iy-1][iz+1] + q[ix+0][iy+1][iz-1] + q[ix+0][iy+1][iz+1] + q[ix+1][iy-1][iz+0] + q[ix+1][iy+0][iz-1] + q[ix+1][iy+0][iz+1] + q[ix+1][iy+1][iz+0])
       + a1*a1*a1*(q[ix-1][iy-1][iz-1] + q[ix-1][iy-1][iz+1] + q[ix-1][iy+1][iz-1] + q[ix-1][iy+1][iz+1] + q[ix+1][iy-1][iz-1] + q[ix+1][iy-1][iz+1] + q[ix+1][iy+1][iz-1] + q[ix+1][iy+1][iz+1])
       + a2*a0*a0*(q[ix-2][iy+0][iz+0] + q[ix+0][iy-2][iz+0] + q[ix+0][iy+0][iz-2] + q[ix+0][iy+0][iz+2] + q[ix+0][iy+2][iz+0] + q[ix+2][iy+0][iz+0])
       + a2*a1*a0*(q[ix-2][iy-1][iz+0] + q[ix-2][iy+0][iz-1] + q[ix-2][iy+0][iz+1] + q[ix-2][iy+1][iz+0] + q[ix-1][iy-2][iz+0] + q[ix-1][iy+0][iz-2] + q[ix-1][iy+0][iz+2] + q[ix-1][iy+2][iz+0] + q[ix+0][iy-2][iz-1] + q[ix+0][iy-2][iz+1] + q[ix+0][iy-1][iz-2] + q[ix+0][iy-1][iz+2] + q[ix+0][iy+1][iz-2] + q[ix+0][iy+1][iz+2] + q[ix+0][iy+2][iz-1] + q[ix+0][iy+2][iz+1] + q[ix+1][iy-2][iz+0] + q[ix+1][iy+0][iz-2] + q[ix+1][iy+0][iz+2] + q[ix+1][iy+2][iz+0] + q[ix+2][iy-1][iz+0] + q[ix+2][iy+0][iz-1] + q[ix+2][iy+0][iz+1] + q[ix+2][iy+1][iz+0])
       + a2*a1*a1*(q[ix-2][iy-1][iz-1] + q[ix-2][iy-1][iz+1] + q[ix-2][iy+1][iz-1] + q[ix-2][iy+1][iz+1] + q[ix-1][iy-2][iz-1] + q[ix-1][iy-2][iz+1] + q[ix-1][iy-1][iz-2] + q[ix-1][iy-1][iz+2] + q[ix-1][iy+1][iz-2] + q[ix-1][iy+1][iz+2] + q[ix-1][iy+2][iz-1] + q[ix-1][iy+2][iz+1] + q[ix+1][iy-2][iz-1] + q[ix+1][iy-2][iz+1] + q[ix+1][iy-1][iz-2] + q[ix+1][iy-1][iz+2] + q[ix+1][iy+1][iz-2] + q[ix+1][iy+1][iz+2] + q[ix+1][iy+2][iz-1] + q[ix+1][iy+2][iz+1] + q[ix+2][iy-1][iz-1] + q[ix+2][iy-1][iz+1] + q[ix+2][iy+1][iz-1] + q[ix+2][iy+1][iz+1])
       + a2*a2*a0*(q[ix-2][iy-2][iz+0] + q[ix-2][iy+0][iz-2] + q[ix-2][iy+0][iz+2] + q[ix-2][iy+2][iz+0] + q[ix+0][iy-2][iz-2] + q[ix+0][iy-2][iz+2] + q[ix+0][iy+2][iz-2] + q[ix+0][iy+2][iz+2] + q[ix+2][iy-2][iz+0] + q[ix+2][iy+0][iz-2] + q[ix+2][iy+0][iz+2] + q[ix+2][iy+2][iz+0])
       + a3*a0*a0*(q[ix-3][iy+0][iz+0] + q[ix+0][iy-3][iz+0] + q[ix+0][iy+0][iz-3] + q[ix+0][iy+0][iz+3] + q[ix+0][iy+3][iz+0] + q[ix+3][iy+0][iz+0])
       + a2*a2*a1*(q[ix-2][iy-2][iz-1] + q[ix-2][iy-2][iz+1] + q[ix-2][iy-1][iz-2] + q[ix-2][iy-1][iz+2] + q[ix-2][iy+1][iz-2] + q[ix-2][iy+1][iz+2] + q[ix-2][iy+2][iz-1] + q[ix-2][iy+2][iz+1] + q[ix-1][iy-2][iz-2] + q[ix-1][iy-2][iz+2] + q[ix-1][iy+2][iz-2] + q[ix-1][iy+2][iz+2] + q[ix+1][iy-2][iz-2] + q[ix+1][iy-2][iz+2] + q[ix+1][iy+2][iz-2] + q[ix+1][iy+2][iz+2] + q[ix+2][iy-2][iz-1] + q[ix+2][iy-2][iz+1] + q[ix+2][iy-1][iz-2] + q[ix+2][iy-1][iz+2] + q[ix+2][iy+1][iz-2] + q[ix+2][iy+1][iz+2] + q[ix+2][iy+2][iz-1] + q[ix+2][iy+2][iz+1])
       + a3*a1*a0*(q[ix-3][iy-1][iz+0] + q[ix-3][iy+0][iz-1] + q[ix-3][iy+0][iz+1] + q[ix-3][iy+1][iz+0] + q[ix-1][iy-3][iz+0] + q[ix-1][iy+0][iz-3] + q[ix-1][iy+0][iz+3] + q[ix-1][iy+3][iz+0] + q[ix+0][iy-3][iz-1] + q[ix+0][iy-3][iz+1] + q[ix+0][iy-1][iz-3] + q[ix+0][iy-1][iz+3] + q[ix+0][iy+1][iz-3] + q[ix+0][iy+1][iz+3] + q[ix+0][iy+3][iz-1] + q[ix+0][iy+3][iz+1] + q[ix+1][iy-3][iz+0] + q[ix+1][iy+0][iz-3] + q[ix+1][iy+0][iz+3] + q[ix+1][iy+3][iz+0] + q[ix+3][iy-1][iz+0] + q[ix+3][iy+0][iz-1] + q[ix+3][iy+0][iz+1] + q[ix+3][iy+1][iz+0])
       + a3*a1*a1*(q[ix-3][iy-1][iz-1] + q[ix-3][iy-1][iz+1] + q[ix-3][iy+1][iz-1] + q[ix-3][iy+1][iz+1] + q[ix-1][iy-3][iz-1] + q[ix-1][iy-3][iz+1] + q[ix-1][iy-1][iz-3] + q[ix-1][iy-1][iz+3] + q[ix-1][iy+1][iz-3] + q[ix-1][iy+1][iz+3] + q[ix-1][iy+3][iz-1] + q[ix-1][iy+3][iz+1] + q[ix+1][iy-3][iz-1] + q[ix+1][iy-3][iz+1] + q[ix+1][iy-1][iz-3] + q[ix+1][iy-1][iz+3] + q[ix+1][iy+1][iz-3] + q[ix+1][iy+1][iz+3] + q[ix+1][iy+3][iz-1] + q[ix+1][iy+3][iz+1] + q[ix+3][iy-1][iz-1] + q[ix+3][iy-1][iz+1] + q[ix+3][iy+1][iz-1] + q[ix+3][iy+1][iz+1])
       + a2*a2*a2*(q[ix-2][iy-2][iz-2] + q[ix-2][iy-2][iz+2] + q[ix-2][iy+2][iz-2] + q[ix-2][iy+2][iz+2] + q[ix+2][iy-2][iz-2] + q[ix+2][iy-2][iz+2] + q[ix+2][iy+2][iz-2] + q[ix+2][iy+2][iz+2])
       + a3*a2*a0*(q[ix-3][iy-2][iz+0] + q[ix-3][iy+0][iz-2] + q[ix-3][iy+0][iz+2] + q[ix-3][iy+2][iz+0] + q[ix-2][iy-3][iz+0] + q[ix-2][iy+0][iz-3] + q[ix-2][iy+0][iz+3] + q[ix-2][iy+3][iz+0] + q[ix+0][iy-3][iz-2] + q[ix+0][iy-3][iz+2] + q[ix+0][iy-2][iz-3] + q[ix+0][iy-2][iz+3] + q[ix+0][iy+2][iz-3] + q[ix+0][iy+2][iz+3] + q[ix+0][iy+3][iz-2] + q[ix+0][iy+3][iz+2] + q[ix+2][iy-3][iz+0] + q[ix+2][iy+0][iz-3] + q[ix+2][iy+0][iz+3] + q[ix+2][iy+3][iz+0] + q[ix+3][iy-2][iz+0] + q[ix+3][iy+0][iz-2] + q[ix+3][iy+0][iz+2] + q[ix+3][iy+2][iz+0])
       + a3*a2*a1*(q[ix-3][iy-2][iz-1] + q[ix-3][iy-2][iz+1] + q[ix-3][iy-1][iz-2] + q[ix-3][iy-1][iz+2] + q[ix-3][iy+1][iz-2] + q[ix-3][iy+1][iz+2] + q[ix-3][iy+2][iz-1] + q[ix-3][iy+2][iz+1] + q[ix-2][iy-3][iz-1] + q[ix-2][iy-3][iz+1] + q[ix-2][iy-1][iz-3] + q[ix-2][iy-1][iz+3] + q[ix-2][iy+1][iz-3] + q[ix-2][iy+1][iz+3] + q[ix-2][iy+3][iz-1] + q[ix-2][iy+3][iz+1] + q[ix-1][iy-3][iz-2] + q[ix-1][iy-3][iz+2] + q[ix-1][iy-2][iz-3] + q[ix-1][iy-2][iz+3] + q[ix-1][iy+2][iz-3] + q[ix-1][iy+2][iz+3] + q[ix-1][iy+3][iz-2] + q[ix-1][iy+3][iz+2] + q[ix+1][iy-3][iz-2] + q[ix+1][iy-3][iz+2] + q[ix+1][iy-2][iz-3] + q[ix+1][iy-2][iz+3] + q[ix+1][iy+2][iz-3] + q[ix+1][iy+2][iz+3] + q[ix+1][iy+3][iz-2] + q[ix+1][iy+3][iz+2] + q[ix+2][iy-3][iz-1] + q[ix+2][iy-3][iz+1] + q[ix+2][iy-1][iz-3] + q[ix+2][iy-1][iz+3] + q[ix+2][iy+1][iz-3] + q[ix+2][iy+1][iz+3] + q[ix+2][iy+3][iz-1] + q[ix+2][iy+3][iz+1] + q[ix+3][iy-2][iz-1] + q[ix+3][iy-2][iz+1] + q[ix+3][iy-1][iz-2] + q[ix+3][iy-1][iz+2] + q[ix+3][iy+1][iz-2] + q[ix+3][iy+1][iz+2] + q[ix+3][iy+2][iz-1] + q[ix+3][iy+2][iz+1])
       + a3*a2*a2*(q[ix-3][iy-2][iz-2] + q[ix-3][iy-2][iz+2] + q[ix-3][iy+2][iz-2] + q[ix-3][iy+2][iz+2] + q[ix-2][iy-3][iz-2] + q[ix-2][iy-3][iz+2] + q[ix-2][iy-2][iz-3] + q[ix-2][iy-2][iz+3] + q[ix-2][iy+2][iz-3] + q[ix-2][iy+2][iz+3] + q[ix-2][iy+3][iz-2] + q[ix-2][iy+3][iz+2] + q[ix+2][iy-3][iz-2] + q[ix+2][iy-3][iz+2] + q[ix+2][iy-2][iz-3] + q[ix+2][iy-2][iz+3] + q[ix+2][iy+2][iz-3] + q[ix+2][iy+2][iz+3] + q[ix+2][iy+3][iz-2] + q[ix+2][iy+3][iz+2] + q[ix+3][iy-2][iz-2] + q[ix+3][iy-2][iz+2] + q[ix+3][iy+2][iz-2] + q[ix+3][iy+2][iz+2])
       + a3*a3*a0*(q[ix-3][iy-3][iz+0] + q[ix-3][iy+0][iz-3] + q[ix-3][iy+0][iz+3] + q[ix-3][iy+3][iz+0] + q[ix+0][iy-3][iz-3] + q[ix+0][iy-3][iz+3] + q[ix+0][iy+3][iz-3] + q[ix+0][iy+3][iz+3] + q[ix+3][iy-3][iz+0] + q[ix+3][iy+0][iz-3] + q[ix+3][iy+0][iz+3] + q[ix+3][iy+3][iz+0])
       + a3*a3*a1*(q[ix-3][iy-3][iz-1] + q[ix-3][iy-3][iz+1] + q[ix-3][iy-1][iz-3] + q[ix-3][iy-1][iz+3] + q[ix-3][iy+1][iz-3] + q[ix-3][iy+1][iz+3] + q[ix-3][iy+3][iz-1] + q[ix-3][iy+3][iz+1] + q[ix-1][iy-3][iz-3] + q[ix-1][iy-3][iz+3] + q[ix-1][iy+3][iz-3] + q[ix-1][iy+3][iz+3] + q[ix+1][iy-3][iz-3] + q[ix+1][iy-3][iz+3] + q[ix+1][iy+3][iz-3] + q[ix+1][iy+3][iz+3] + q[ix+3][iy-3][iz-1] + q[ix+3][iy-3][iz+1] + q[ix+3][iy-1][iz-3] + q[ix+3][iy-1][iz+3] + q[ix+3][iy+1][iz-3] + q[ix+3][iy+1][iz+3] + q[ix+3][iy+3][iz-1] + q[ix+3][iy+3][iz+1])
       + a3*a3*a2*(q[ix-3][iy-3][iz-2] + q[ix-3][iy-3][iz+2] + q[ix-3][iy-2][iz-3] + q[ix-3][iy-2][iz+3] + q[ix-3][iy+2][iz-3] + q[ix-3][iy+2][iz+3] + q[ix-3][iy+3][iz-2] + q[ix-3][iy+3][iz+2] + q[ix-2][iy-3][iz-3] + q[ix-2][iy-3][iz+3] + q[ix-2][iy+3][iz-3] + q[ix-2][iy+3][iz+3] + q[ix+2][iy-3][iz-3] + q[ix+2][iy-3][iz+3] + q[ix+2][iy+3][iz-3] + q[ix+2][iy+3][iz+3] + q[ix+3][iy-3][iz-2] + q[ix+3][iy-3][iz+2] + q[ix+3][iy-2][iz-3] + q[ix+3][iy-2][iz+3] + q[ix+3][iy+2][iz-3] + q[ix+3][iy+2][iz+3] + q[ix+3][iy+3][iz-2] + q[ix+3][iy+3][iz+2])
       + a3*a3*a3*(q[ix-3][iy-3][iz-3] + q[ix-3][iy-3][iz+3] + q[ix-3][iy+3][iz-3] + q[ix-3][iy+3][iz+3] + q[ix+3][iy-3][iz-3] + q[ix+3][iy-3][iz+3] + q[ix+3][iy+3][iz-3] + q[ix+3][iy+3][iz+3]);
}

struct tmp_buf {
  double rp_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double up_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double vp_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double wp_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double pp_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];

  double rh_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double uh_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double vh_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double wh_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
  double ph_tmp_buf[MX*NX+6][MY*NY+6][MZ*NZ+6];
};

void lowpass_filter(navi &n, floors<double> &s) {
/*

  static auto tmp = unique_ptr<tmp_buf>(new tmp_buf);

  // 床の準備 (本来はMPI通信が必要)
  for(int ix = 0; ix < MX*NX+2*Ms; ix++) {
    int x = (ix-Ms+LX)%LX + 2*NX;
    for(int iy = 0; iy < MY*NY+2*Ms; iy++) {
      int y = (iy-Ms+LY)%LY + 2*NY;
      for(int iz = 0; iz < MZ*NZ+2*Ms; iz++) {
        int z = (iz-Ms+LZ)%LZ + 2*NZ;
        tmp->rp_tmp_buf[ix][iy][iz] = s.p[x][y][z].r;
        tmp->up_tmp_buf[ix][iy][iz] = s.p[x][y][z].u;
        tmp->vp_tmp_buf[ix][iy][iz] = s.p[x][y][z].v;
        tmp->wp_tmp_buf[ix][iy][iz] = s.p[x][y][z].w;
        tmp->pp_tmp_buf[ix][iy][iz] = s.p[x][y][z].p;

        tmp->rh_tmp_buf[ix][iy][iz] = s.h[x][y][z].r;
        tmp->uh_tmp_buf[ix][iy][iz] = s.h[x][y][z].u;
        tmp->vh_tmp_buf[ix][iy][iz] = s.h[x][y][z].v;
        tmp->wh_tmp_buf[ix][iy][iz] = s.h[x][y][z].w;
        tmp->ph_tmp_buf[ix][iy][iz] = s.h[x][y][z].p;
      }
    }
  }

  // フィルターの適用
  for(int ix = Ms; ix < MX*NX+Ms; ix++) {
    for(int iy = Ms; iy < MY*NY+Ms; iy++) {
      for(int iz = Ms; iz < MZ*NZ+Ms; iz++) {
        tmp->rp_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->rp_tmp_buf);
        tmp->up_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->up_tmp_buf);
        tmp->vp_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->vp_tmp_buf);
        tmp->wp_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->wp_tmp_buf);
        tmp->pp_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->pp_tmp_buf);

        tmp->rh_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->rh_tmp_buf);
        tmp->uh_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->uh_tmp_buf);
        tmp->vh_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->vh_tmp_buf);
        tmp->wh_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->wh_tmp_buf);
        tmp->ph_tmp_buf[ix-Ms][iy-Ms][iz-Ms] = filter(ix,iy,iz,tmp->ph_tmp_buf);
      }
    }
  }

  // 床の書きこみ
  for(int ix = 0; ix < MX*NX; ix++) {
    int x = ix + 2*NX;
    for(int iy = 0; iy < MY*NY; iy++) {
      int y = iy + 2*NY;
      for(int iz = 0; iz < MZ*NZ; iz++) {
        int z = iz + 2*NZ;
        s.p[x][y][z].r = tmp->rp_tmp_buf[ix][iy][iz];
        s.p[x][y][z].u = tmp->up_tmp_buf[ix][iy][iz];
        s.p[x][y][z].v = tmp->vp_tmp_buf[ix][iy][iz];
        s.p[x][y][z].w = tmp->wp_tmp_buf[ix][iy][iz];
        s.p[x][y][z].p = tmp->pp_tmp_buf[ix][iy][iz];

        s.h[x][y][z].r = tmp->rh_tmp_buf[ix][iy][iz];
        s.h[x][y][z].u = tmp->uh_tmp_buf[ix][iy][iz];
        s.h[x][y][z].v = tmp->vh_tmp_buf[ix][iy][iz];
        s.h[x][y][z].w = tmp->wh_tmp_buf[ix][iy][iz];
        s.h[x][y][z].p = tmp->ph_tmp_buf[ix][iy][iz];
      }
    }
  }
  */
}
