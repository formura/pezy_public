#include <math.h>
#include <stdio.h>
#include <memory>
#include "config.h"

using namespace std;

typedef struct {
  // 床の最小座標
  int x_origin;
  int y_origin;
  int z_origin;

  double rp_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double up_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double vp_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double wp_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double pp_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double rp_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double up_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double vp_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double wp_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double pp_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double rp_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double up_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double vp_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double wp_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double pp_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];

  double rh_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double uh_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double vh_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double wh_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double ph_wall_x[NT][2*Ns][NY+2*Ns][NZ+2*Ns];
  double rh_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double uh_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double vh_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double wh_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double ph_wall_y[NT][NX+2*Ns][2*Ns][NZ+2*Ns];
  double rh_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double uh_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double vh_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double wh_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];
  double ph_wall_z[NT][NX+2*Ns][NY+2*Ns][2*Ns];

  double rp_res[NX][NY][NZ];
  double up_res[NX][NY][NZ];
  double vp_res[NX][NY][NZ];
  double wp_res[NX][NY][NZ];
  double pp_res[NX][NY][NZ];
  double rh_res[NX][NY][NZ];
  double uh_res[NX][NY][NZ];
  double vh_res[NX][NY][NZ];
  double wh_res[NX][NY][NZ];
  double ph_res[NX][NY][NZ];
} param;

double d_x(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j+0][k+0] -8 * q[i-1][j+0][k+0] +8 * q[i+1][j+0][k+0] -1 * q[i+2][j+0][k+0])/(12*h);}
double d_y(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j-2][k+0] -8 * q[i+0][j-1][k+0] +8 * q[i+0][j+1][k+0] -1 * q[i+0][j+2][k+0])/(12*h);}
double d_z(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j+0][k-2] -8 * q[i+0][j+0][k-1] +8 * q[i+0][j+0][k+1] -1 * q[i+0][j+0][k+2])/(12*h);}
double d_xx(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k+0] +16 * q[i-1][j+0][k+0] -30 * q[i+0][j+0][k+0] +16 * q[i+1][j+0][k+0] -1 * q[i+2][j+0][k+0])/(12*h*h);}
double d_xy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j-2][k+0] -8 * q[i-2][j-1][k+0] +8 * q[i-2][j+1][k+0] -1 * q[i-2][j+2][k+0] -8 * q[i-1][j-2][k+0] +64 * q[i-1][j-1][k+0] -64 * q[i-1][j+1][k+0] +8 * q[i-1][j+2][k+0] +8 * q[i+1][j-2][k+0] -64 * q[i+1][j-1][k+0] +64 * q[i+1][j+1][k+0] -8 * q[i+1][j+2][k+0] -1 * q[i+2][j-2][k+0] +8 * q[i+2][j-1][k+0] -8 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h);}
double d_yy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k+0] +16 * q[i+0][j-1][k+0] -30 * q[i+0][j+0][k+0] +16 * q[i+0][j+1][k+0] -1 * q[i+0][j+2][k+0])/(12*h*h);}
double d_xz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j+0][k-2] -8 * q[i-2][j+0][k-1] +8 * q[i-2][j+0][k+1] -1 * q[i-2][j+0][k+2] -8 * q[i-1][j+0][k-2] +64 * q[i-1][j+0][k-1] -64 * q[i-1][j+0][k+1] +8 * q[i-1][j+0][k+2] +8 * q[i+1][j+0][k-2] -64 * q[i+1][j+0][k-1] +64 * q[i+1][j+0][k+1] -8 * q[i+1][j+0][k+2] -1 * q[i+2][j+0][k-2] +8 * q[i+2][j+0][k-1] -8 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h);}
double d_yz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j-2][k-2] -8 * q[i+0][j-2][k-1] +8 * q[i+0][j-2][k+1] -1 * q[i+0][j-2][k+2] -8 * q[i+0][j-1][k-2] +64 * q[i+0][j-1][k-1] -64 * q[i+0][j-1][k+1] +8 * q[i+0][j-1][k+2] +8 * q[i+0][j+1][k-2] -64 * q[i+0][j+1][k-1] +64 * q[i+0][j+1][k+1] -8 * q[i+0][j+1][k+2] -1 * q[i+0][j+2][k-2] +8 * q[i+0][j+2][k-1] -8 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h);}
double d_zz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j+0][k-2] +16 * q[i+0][j+0][k-1] -30 * q[i+0][j+0][k+0] +16 * q[i+0][j+0][k+1] -1 * q[i+0][j+0][k+2])/(12*h*h);}
double d_xxx(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k+0] +2 * q[i-1][j+0][k+0] -2 * q[i+1][j+0][k+0] +1 * q[i+2][j+0][k+0])/(2*h*h*h);}
double d_xxy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +8 * q[i-2][j-1][k+0] -8 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +16 * q[i-1][j-2][k+0] -128 * q[i-1][j-1][k+0] +128 * q[i-1][j+1][k+0] -16 * q[i-1][j+2][k+0] -30 * q[i+0][j-2][k+0] +240 * q[i+0][j-1][k+0] -240 * q[i+0][j+1][k+0] +30 * q[i+0][j+2][k+0] +16 * q[i+1][j-2][k+0] -128 * q[i+1][j-1][k+0] +128 * q[i+1][j+1][k+0] -16 * q[i+1][j+2][k+0] -1 * q[i+2][j-2][k+0] +8 * q[i+2][j-1][k+0] -8 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h);}
double d_xyy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +16 * q[i-2][j-1][k+0] -30 * q[i-2][j+0][k+0] +16 * q[i-2][j+1][k+0] -1 * q[i-2][j+2][k+0] +8 * q[i-1][j-2][k+0] -128 * q[i-1][j-1][k+0] +240 * q[i-1][j+0][k+0] -128 * q[i-1][j+1][k+0] +8 * q[i-1][j+2][k+0] -8 * q[i+1][j-2][k+0] +128 * q[i+1][j-1][k+0] -240 * q[i+1][j+0][k+0] +128 * q[i+1][j+1][k+0] -8 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -16 * q[i+2][j-1][k+0] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h);}
double d_yyy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k+0] +2 * q[i+0][j-1][k+0] -2 * q[i+0][j+1][k+0] +1 * q[i+0][j+2][k+0])/(2*h*h*h);}
double d_xxz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +8 * q[i-2][j+0][k-1] -8 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +16 * q[i-1][j+0][k-2] -128 * q[i-1][j+0][k-1] +128 * q[i-1][j+0][k+1] -16 * q[i-1][j+0][k+2] -30 * q[i+0][j+0][k-2] +240 * q[i+0][j+0][k-1] -240 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] +16 * q[i+1][j+0][k-2] -128 * q[i+1][j+0][k-1] +128 * q[i+1][j+0][k+1] -16 * q[i+1][j+0][k+2] -1 * q[i+2][j+0][k-2] +8 * q[i+2][j+0][k-1] -8 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h);}
double d_xyz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j-2][k-2] -8 * q[i-2][j-2][k-1] +8 * q[i-2][j-2][k+1] -1 * q[i-2][j-2][k+2] -8 * q[i-2][j-1][k-2] +64 * q[i-2][j-1][k-1] -64 * q[i-2][j-1][k+1] +8 * q[i-2][j-1][k+2] +8 * q[i-2][j+1][k-2] -64 * q[i-2][j+1][k-1] +64 * q[i-2][j+1][k+1] -8 * q[i-2][j+1][k+2] -1 * q[i-2][j+2][k-2] +8 * q[i-2][j+2][k-1] -8 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] -8 * q[i-1][j-2][k-2] +64 * q[i-1][j-2][k-1] -64 * q[i-1][j-2][k+1] +8 * q[i-1][j-2][k+2] +64 * q[i-1][j-1][k-2] -512 * q[i-1][j-1][k-1] +512 * q[i-1][j-1][k+1] -64 * q[i-1][j-1][k+2] -64 * q[i-1][j+1][k-2] +512 * q[i-1][j+1][k-1] -512 * q[i-1][j+1][k+1] +64 * q[i-1][j+1][k+2] +8 * q[i-1][j+2][k-2] -64 * q[i-1][j+2][k-1] +64 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] +8 * q[i+1][j-2][k-2] -64 * q[i+1][j-2][k-1] +64 * q[i+1][j-2][k+1] -8 * q[i+1][j-2][k+2] -64 * q[i+1][j-1][k-2] +512 * q[i+1][j-1][k-1] -512 * q[i+1][j-1][k+1] +64 * q[i+1][j-1][k+2] +64 * q[i+1][j+1][k-2] -512 * q[i+1][j+1][k-1] +512 * q[i+1][j+1][k+1] -64 * q[i+1][j+1][k+2] -8 * q[i+1][j+2][k-2] +64 * q[i+1][j+2][k-1] -64 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] -1 * q[i+2][j-2][k-2] +8 * q[i+2][j-2][k-1] -8 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] +8 * q[i+2][j-1][k-2] -64 * q[i+2][j-1][k-1] +64 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] -8 * q[i+2][j+1][k-2] +64 * q[i+2][j+1][k-1] -64 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h);}
double d_yyz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +8 * q[i+0][j-2][k-1] -8 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +16 * q[i+0][j-1][k-2] -128 * q[i+0][j-1][k-1] +128 * q[i+0][j-1][k+1] -16 * q[i+0][j-1][k+2] -30 * q[i+0][j+0][k-2] +240 * q[i+0][j+0][k-1] -240 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] +16 * q[i+0][j+1][k-2] -128 * q[i+0][j+1][k-1] +128 * q[i+0][j+1][k+1] -16 * q[i+0][j+1][k+2] -1 * q[i+0][j+2][k-2] +8 * q[i+0][j+2][k-1] -8 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h);}
double d_xzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +16 * q[i-2][j+0][k-1] -30 * q[i-2][j+0][k+0] +16 * q[i-2][j+0][k+1] -1 * q[i-2][j+0][k+2] +8 * q[i-1][j+0][k-2] -128 * q[i-1][j+0][k-1] +240 * q[i-1][j+0][k+0] -128 * q[i-1][j+0][k+1] +8 * q[i-1][j+0][k+2] -8 * q[i+1][j+0][k-2] +128 * q[i+1][j+0][k-1] -240 * q[i+1][j+0][k+0] +128 * q[i+1][j+0][k+1] -8 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -16 * q[i+2][j+0][k-1] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h);}
double d_yzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +16 * q[i+0][j-2][k-1] -30 * q[i+0][j-2][k+0] +16 * q[i+0][j-2][k+1] -1 * q[i+0][j-2][k+2] +8 * q[i+0][j-1][k-2] -128 * q[i+0][j-1][k-1] +240 * q[i+0][j-1][k+0] -128 * q[i+0][j-1][k+1] +8 * q[i+0][j-1][k+2] -8 * q[i+0][j+1][k-2] +128 * q[i+0][j+1][k-1] -240 * q[i+0][j+1][k+0] +128 * q[i+0][j+1][k+1] -8 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -16 * q[i+0][j+2][k-1] +30 * q[i+0][j+2][k+0] -16 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h);}
double d_zzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j+0][k-2] +2 * q[i+0][j+0][k-1] -2 * q[i+0][j+0][k+1] +1 * q[i+0][j+0][k+2])/(2*h*h*h);}
double d_xxxx(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j+0][k+0] -4 * q[i-1][j+0][k+0] +6 * q[i+0][j+0][k+0] -4 * q[i+1][j+0][k+0] +1 * q[i+2][j+0][k+0])/(1*h*h*h*h);}
double d_xxxy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +8 * q[i-2][j-1][k+0] -8 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +2 * q[i-1][j-2][k+0] -16 * q[i-1][j-1][k+0] +16 * q[i-1][j+1][k+0] -2 * q[i-1][j+2][k+0] -2 * q[i+1][j-2][k+0] +16 * q[i+1][j-1][k+0] -16 * q[i+1][j+1][k+0] +2 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -8 * q[i+2][j-1][k+0] +8 * q[i+2][j+1][k+0] -1 * q[i+2][j+2][k+0])/(24*h*h*h*h);}
double d_xxyy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j-2][k+0] -16 * q[i-2][j-1][k+0] +30 * q[i-2][j+0][k+0] -16 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] -16 * q[i-1][j-2][k+0] +256 * q[i-1][j-1][k+0] -480 * q[i-1][j+0][k+0] +256 * q[i-1][j+1][k+0] -16 * q[i-1][j+2][k+0] +30 * q[i+0][j-2][k+0] -480 * q[i+0][j-1][k+0] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+1][k+0] +30 * q[i+0][j+2][k+0] -16 * q[i+1][j-2][k+0] +256 * q[i+1][j-1][k+0] -480 * q[i+1][j+0][k+0] +256 * q[i+1][j+1][k+0] -16 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -16 * q[i+2][j-1][k+0] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+1][k+0] +1 * q[i+2][j+2][k+0])/(144*h*h*h*h);}
double d_xyyy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k+0] +2 * q[i-2][j-1][k+0] -2 * q[i-2][j+1][k+0] +1 * q[i-2][j+2][k+0] +8 * q[i-1][j-2][k+0] -16 * q[i-1][j-1][k+0] +16 * q[i-1][j+1][k+0] -8 * q[i-1][j+2][k+0] -8 * q[i+1][j-2][k+0] +16 * q[i+1][j-1][k+0] -16 * q[i+1][j+1][k+0] +8 * q[i+1][j+2][k+0] +1 * q[i+2][j-2][k+0] -2 * q[i+2][j-1][k+0] +2 * q[i+2][j+1][k+0] -1 * q[i+2][j+2][k+0])/(24*h*h*h*h);}
double d_yyyy(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j-2][k+0] -4 * q[i+0][j-1][k+0] +6 * q[i+0][j+0][k+0] -4 * q[i+0][j+1][k+0] +1 * q[i+0][j+2][k+0])/(1*h*h*h*h);}
double d_xxxz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +8 * q[i-2][j+0][k-1] -8 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +2 * q[i-1][j+0][k-2] -16 * q[i-1][j+0][k-1] +16 * q[i-1][j+0][k+1] -2 * q[i-1][j+0][k+2] -2 * q[i+1][j+0][k-2] +16 * q[i+1][j+0][k-1] -16 * q[i+1][j+0][k+1] +2 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -8 * q[i+2][j+0][k-1] +8 * q[i+2][j+0][k+1] -1 * q[i+2][j+0][k+2])/(24*h*h*h*h);}
double d_xxyz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +8 * q[i-2][j-2][k-1] -8 * q[i-2][j-2][k+1] +1 * q[i-2][j-2][k+2] +8 * q[i-2][j-1][k-2] -64 * q[i-2][j-1][k-1] +64 * q[i-2][j-1][k+1] -8 * q[i-2][j-1][k+2] -8 * q[i-2][j+1][k-2] +64 * q[i-2][j+1][k-1] -64 * q[i-2][j+1][k+1] +8 * q[i-2][j+1][k+2] +1 * q[i-2][j+2][k-2] -8 * q[i-2][j+2][k-1] +8 * q[i-2][j+2][k+1] -1 * q[i-2][j+2][k+2] +16 * q[i-1][j-2][k-2] -128 * q[i-1][j-2][k-1] +128 * q[i-1][j-2][k+1] -16 * q[i-1][j-2][k+2] -128 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1024 * q[i-1][j-1][k+1] +128 * q[i-1][j-1][k+2] +128 * q[i-1][j+1][k-2] -1024 * q[i-1][j+1][k-1] +1024 * q[i-1][j+1][k+1] -128 * q[i-1][j+1][k+2] -16 * q[i-1][j+2][k-2] +128 * q[i-1][j+2][k-1] -128 * q[i-1][j+2][k+1] +16 * q[i-1][j+2][k+2] -30 * q[i+0][j-2][k-2] +240 * q[i+0][j-2][k-1] -240 * q[i+0][j-2][k+1] +30 * q[i+0][j-2][k+2] +240 * q[i+0][j-1][k-2] -1920 * q[i+0][j-1][k-1] +1920 * q[i+0][j-1][k+1] -240 * q[i+0][j-1][k+2] -240 * q[i+0][j+1][k-2] +1920 * q[i+0][j+1][k-1] -1920 * q[i+0][j+1][k+1] +240 * q[i+0][j+1][k+2] +30 * q[i+0][j+2][k-2] -240 * q[i+0][j+2][k-1] +240 * q[i+0][j+2][k+1] -30 * q[i+0][j+2][k+2] +16 * q[i+1][j-2][k-2] -128 * q[i+1][j-2][k-1] +128 * q[i+1][j-2][k+1] -16 * q[i+1][j-2][k+2] -128 * q[i+1][j-1][k-2] +1024 * q[i+1][j-1][k-1] -1024 * q[i+1][j-1][k+1] +128 * q[i+1][j-1][k+2] +128 * q[i+1][j+1][k-2] -1024 * q[i+1][j+1][k-1] +1024 * q[i+1][j+1][k+1] -128 * q[i+1][j+1][k+2] -16 * q[i+1][j+2][k-2] +128 * q[i+1][j+2][k-1] -128 * q[i+1][j+2][k+1] +16 * q[i+1][j+2][k+2] -1 * q[i+2][j-2][k-2] +8 * q[i+2][j-2][k-1] -8 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] +8 * q[i+2][j-1][k-2] -64 * q[i+2][j-1][k-1] +64 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] -8 * q[i+2][j+1][k-2] +64 * q[i+2][j+1][k-1] -64 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_xyyz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +8 * q[i-2][j-2][k-1] -8 * q[i-2][j-2][k+1] +1 * q[i-2][j-2][k+2] +16 * q[i-2][j-1][k-2] -128 * q[i-2][j-1][k-1] +128 * q[i-2][j-1][k+1] -16 * q[i-2][j-1][k+2] -30 * q[i-2][j+0][k-2] +240 * q[i-2][j+0][k-1] -240 * q[i-2][j+0][k+1] +30 * q[i-2][j+0][k+2] +16 * q[i-2][j+1][k-2] -128 * q[i-2][j+1][k-1] +128 * q[i-2][j+1][k+1] -16 * q[i-2][j+1][k+2] -1 * q[i-2][j+2][k-2] +8 * q[i-2][j+2][k-1] -8 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] +8 * q[i-1][j-2][k-2] -64 * q[i-1][j-2][k-1] +64 * q[i-1][j-2][k+1] -8 * q[i-1][j-2][k+2] -128 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1024 * q[i-1][j-1][k+1] +128 * q[i-1][j-1][k+2] +240 * q[i-1][j+0][k-2] -1920 * q[i-1][j+0][k-1] +1920 * q[i-1][j+0][k+1] -240 * q[i-1][j+0][k+2] -128 * q[i-1][j+1][k-2] +1024 * q[i-1][j+1][k-1] -1024 * q[i-1][j+1][k+1] +128 * q[i-1][j+1][k+2] +8 * q[i-1][j+2][k-2] -64 * q[i-1][j+2][k-1] +64 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] -8 * q[i+1][j-2][k-2] +64 * q[i+1][j-2][k-1] -64 * q[i+1][j-2][k+1] +8 * q[i+1][j-2][k+2] +128 * q[i+1][j-1][k-2] -1024 * q[i+1][j-1][k-1] +1024 * q[i+1][j-1][k+1] -128 * q[i+1][j-1][k+2] -240 * q[i+1][j+0][k-2] +1920 * q[i+1][j+0][k-1] -1920 * q[i+1][j+0][k+1] +240 * q[i+1][j+0][k+2] +128 * q[i+1][j+1][k-2] -1024 * q[i+1][j+1][k-1] +1024 * q[i+1][j+1][k+1] -128 * q[i+1][j+1][k+2] -8 * q[i+1][j+2][k-2] +64 * q[i+1][j+2][k-1] -64 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] +1 * q[i+2][j-2][k-2] -8 * q[i+2][j-2][k-1] +8 * q[i+2][j-2][k+1] -1 * q[i+2][j-2][k+2] -16 * q[i+2][j-1][k-2] +128 * q[i+2][j-1][k-1] -128 * q[i+2][j-1][k+1] +16 * q[i+2][j-1][k+2] +30 * q[i+2][j+0][k-2] -240 * q[i+2][j+0][k-1] +240 * q[i+2][j+0][k+1] -30 * q[i+2][j+0][k+2] -16 * q[i+2][j+1][k-2] +128 * q[i+2][j+1][k-1] -128 * q[i+2][j+1][k+1] +16 * q[i+2][j+1][k+2] +1 * q[i+2][j+2][k-2] -8 * q[i+2][j+2][k-1] +8 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_yyyz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +8 * q[i+0][j-2][k-1] -8 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +2 * q[i+0][j-1][k-2] -16 * q[i+0][j-1][k-1] +16 * q[i+0][j-1][k+1] -2 * q[i+0][j-1][k+2] -2 * q[i+0][j+1][k-2] +16 * q[i+0][j+1][k-1] -16 * q[i+0][j+1][k+1] +2 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -8 * q[i+0][j+2][k-1] +8 * q[i+0][j+2][k+1] -1 * q[i+0][j+2][k+2])/(24*h*h*h*h);}
double d_xxzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i-2][j+0][k-2] -16 * q[i-2][j+0][k-1] +30 * q[i-2][j+0][k+0] -16 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] -16 * q[i-1][j+0][k-2] +256 * q[i-1][j+0][k-1] -480 * q[i-1][j+0][k+0] +256 * q[i-1][j+0][k+1] -16 * q[i-1][j+0][k+2] +30 * q[i+0][j+0][k-2] -480 * q[i+0][j+0][k-1] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] -16 * q[i+1][j+0][k-2] +256 * q[i+1][j+0][k-1] -480 * q[i+1][j+0][k+0] +256 * q[i+1][j+0][k+1] -16 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -16 * q[i+2][j+0][k-1] +30 * q[i+2][j+0][k+0] -16 * q[i+2][j+0][k+1] +1 * q[i+2][j+0][k+2])/(144*h*h*h*h);}
double d_xyzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j-2][k-2] +16 * q[i-2][j-2][k-1] -30 * q[i-2][j-2][k+0] +16 * q[i-2][j-2][k+1] -1 * q[i-2][j-2][k+2] +8 * q[i-2][j-1][k-2] -128 * q[i-2][j-1][k-1] +240 * q[i-2][j-1][k+0] -128 * q[i-2][j-1][k+1] +8 * q[i-2][j-1][k+2] -8 * q[i-2][j+1][k-2] +128 * q[i-2][j+1][k-1] -240 * q[i-2][j+1][k+0] +128 * q[i-2][j+1][k+1] -8 * q[i-2][j+1][k+2] +1 * q[i-2][j+2][k-2] -16 * q[i-2][j+2][k-1] +30 * q[i-2][j+2][k+0] -16 * q[i-2][j+2][k+1] +1 * q[i-2][j+2][k+2] +8 * q[i-1][j-2][k-2] -128 * q[i-1][j-2][k-1] +240 * q[i-1][j-2][k+0] -128 * q[i-1][j-2][k+1] +8 * q[i-1][j-2][k+2] -64 * q[i-1][j-1][k-2] +1024 * q[i-1][j-1][k-1] -1920 * q[i-1][j-1][k+0] +1024 * q[i-1][j-1][k+1] -64 * q[i-1][j-1][k+2] +64 * q[i-1][j+1][k-2] -1024 * q[i-1][j+1][k-1] +1920 * q[i-1][j+1][k+0] -1024 * q[i-1][j+1][k+1] +64 * q[i-1][j+1][k+2] -8 * q[i-1][j+2][k-2] +128 * q[i-1][j+2][k-1] -240 * q[i-1][j+2][k+0] +128 * q[i-1][j+2][k+1] -8 * q[i-1][j+2][k+2] -8 * q[i+1][j-2][k-2] +128 * q[i+1][j-2][k-1] -240 * q[i+1][j-2][k+0] +128 * q[i+1][j-2][k+1] -8 * q[i+1][j-2][k+2] +64 * q[i+1][j-1][k-2] -1024 * q[i+1][j-1][k-1] +1920 * q[i+1][j-1][k+0] -1024 * q[i+1][j-1][k+1] +64 * q[i+1][j-1][k+2] -64 * q[i+1][j+1][k-2] +1024 * q[i+1][j+1][k-1] -1920 * q[i+1][j+1][k+0] +1024 * q[i+1][j+1][k+1] -64 * q[i+1][j+1][k+2] +8 * q[i+1][j+2][k-2] -128 * q[i+1][j+2][k-1] +240 * q[i+1][j+2][k+0] -128 * q[i+1][j+2][k+1] +8 * q[i+1][j+2][k+2] +1 * q[i+2][j-2][k-2] -16 * q[i+2][j-2][k-1] +30 * q[i+2][j-2][k+0] -16 * q[i+2][j-2][k+1] +1 * q[i+2][j-2][k+2] -8 * q[i+2][j-1][k-2] +128 * q[i+2][j-1][k-1] -240 * q[i+2][j-1][k+0] +128 * q[i+2][j-1][k+1] -8 * q[i+2][j-1][k+2] +8 * q[i+2][j+1][k-2] -128 * q[i+2][j+1][k-1] +240 * q[i+2][j+1][k+0] -128 * q[i+2][j+1][k+1] +8 * q[i+2][j+1][k+2] -1 * q[i+2][j+2][k-2] +16 * q[i+2][j+2][k-1] -30 * q[i+2][j+2][k+0] +16 * q[i+2][j+2][k+1] -1 * q[i+2][j+2][k+2])/(1728*h*h*h*h);}
double d_yyzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j-2][k-2] -16 * q[i+0][j-2][k-1] +30 * q[i+0][j-2][k+0] -16 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] -16 * q[i+0][j-1][k-2] +256 * q[i+0][j-1][k-1] -480 * q[i+0][j-1][k+0] +256 * q[i+0][j-1][k+1] -16 * q[i+0][j-1][k+2] +30 * q[i+0][j+0][k-2] -480 * q[i+0][j+0][k-1] +900 * q[i+0][j+0][k+0] -480 * q[i+0][j+0][k+1] +30 * q[i+0][j+0][k+2] -16 * q[i+0][j+1][k-2] +256 * q[i+0][j+1][k-1] -480 * q[i+0][j+1][k+0] +256 * q[i+0][j+1][k+1] -16 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -16 * q[i+0][j+2][k-1] +30 * q[i+0][j+2][k+0] -16 * q[i+0][j+2][k+1] +1 * q[i+0][j+2][k+2])/(144*h*h*h*h);}
double d_xzzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i-2][j+0][k-2] +2 * q[i-2][j+0][k-1] -2 * q[i-2][j+0][k+1] +1 * q[i-2][j+0][k+2] +8 * q[i-1][j+0][k-2] -16 * q[i-1][j+0][k-1] +16 * q[i-1][j+0][k+1] -8 * q[i-1][j+0][k+2] -8 * q[i+1][j+0][k-2] +16 * q[i+1][j+0][k-1] -16 * q[i+1][j+0][k+1] +8 * q[i+1][j+0][k+2] +1 * q[i+2][j+0][k-2] -2 * q[i+2][j+0][k-1] +2 * q[i+2][j+0][k+1] -1 * q[i+2][j+0][k+2])/(24*h*h*h*h);}
double d_yzzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (-1 * q[i+0][j-2][k-2] +2 * q[i+0][j-2][k-1] -2 * q[i+0][j-2][k+1] +1 * q[i+0][j-2][k+2] +8 * q[i+0][j-1][k-2] -16 * q[i+0][j-1][k-1] +16 * q[i+0][j-1][k+1] -8 * q[i+0][j-1][k+2] -8 * q[i+0][j+1][k-2] +16 * q[i+0][j+1][k-1] -16 * q[i+0][j+1][k+1] +8 * q[i+0][j+1][k+2] +1 * q[i+0][j+2][k-2] -2 * q[i+0][j+2][k-1] +2 * q[i+0][j+2][k+1] -1 * q[i+0][j+2][k+2])/(24*h*h*h*h);}
double d_zzzz(int i,int j,int k,double q[NX+2*Ns][NY+2*Ns][NZ+2*Ns]) { return (+1 * q[i+0][j+0][k-2] -4 * q[i+0][j+0][k-1] +6 * q[i+0][j+0][k+0] -4 * q[i+0][j+0][k+1] +1 * q[i+0][j+0][k+2])/(1*h*h*h*h);}

// 1ステップ更新
void step(double rp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double up_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double vp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double wp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double pp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double rh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double uh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double vh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double wh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns], double ph_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns])
{
  for(int ix = Ns; ix < NX + Ns; ix++) {
    for(int iy = Ns; iy < NY + Ns; iy++) {
      for(int iz = Ns; iz < NZ + Ns; iz++) {
        double r = rp_buf[ix][iy][iz];
        double u = up_buf[ix][iy][iz];
        double v = vp_buf[ix][iy][iz];
        double w = wp_buf[ix][iy][iz];
        double p = pp_buf[ix][iy][iz];

        double r_x = d_x(ix,iy,iz,rp_buf);
        double u_x = d_x(ix,iy,iz,up_buf);
        double v_x = d_x(ix,iy,iz,vp_buf);
        double w_x = d_x(ix,iy,iz,wp_buf);
        double p_x = d_x(ix,iy,iz,pp_buf);
        double r_y = d_y(ix,iy,iz,rp_buf);
        double u_y = d_y(ix,iy,iz,up_buf);
        double v_y = d_y(ix,iy,iz,vp_buf);
        double w_y = d_y(ix,iy,iz,wp_buf);
        double p_y = d_y(ix,iy,iz,pp_buf);
        double r_z = d_z(ix,iy,iz,rp_buf);
        double u_z = d_z(ix,iy,iz,up_buf);
        double v_z = d_z(ix,iy,iz,vp_buf);
        double w_z = d_z(ix,iy,iz,wp_buf);
        double p_z = d_z(ix,iy,iz,pp_buf);

        double r_xx = d_xx(ix,iy,iz,rp_buf);
        double u_xx = d_xx(ix,iy,iz,up_buf);
        double v_xx = d_xx(ix,iy,iz,vp_buf);
        double w_xx = d_xx(ix,iy,iz,wp_buf);
        double p_xx = d_xx(ix,iy,iz,pp_buf);
        double r_yy = d_yy(ix,iy,iz,rp_buf);
        double u_yy = d_yy(ix,iy,iz,up_buf);
        double v_yy = d_yy(ix,iy,iz,vp_buf);
        double w_yy = d_yy(ix,iy,iz,wp_buf);
        double p_yy = d_yy(ix,iy,iz,pp_buf);
        double r_zz = d_zz(ix,iy,iz,rp_buf);
        double u_zz = d_zz(ix,iy,iz,up_buf);
        double v_zz = d_zz(ix,iy,iz,vp_buf);
        double w_zz = d_zz(ix,iy,iz,wp_buf);
        double p_zz = d_zz(ix,iy,iz,pp_buf);
        double r_xy = d_xy(ix,iy,iz,rp_buf);
        double u_xy = d_xy(ix,iy,iz,up_buf);
        double v_xy = d_xy(ix,iy,iz,vp_buf);
        double w_xy = d_xy(ix,iy,iz,wp_buf);
        double p_xy = d_xy(ix,iy,iz,pp_buf);
        double r_yz = d_yz(ix,iy,iz,rp_buf);
        double u_yz = d_yz(ix,iy,iz,up_buf);
        double v_yz = d_yz(ix,iy,iz,vp_buf);
        double w_yz = d_yz(ix,iy,iz,wp_buf);
        double p_yz = d_yz(ix,iy,iz,pp_buf);
        double r_xz = d_xz(ix,iy,iz,rp_buf);
        double u_xz = d_xz(ix,iy,iz,up_buf);
        double v_xz = d_xz(ix,iy,iz,vp_buf);
        double w_xz = d_xz(ix,iy,iz,wp_buf);
        double p_xz = d_xz(ix,iy,iz,pp_buf);

        double u_xxx = d_xxx(ix,iy,iz,up_buf);
        double v_xxx = d_xxx(ix,iy,iz,vp_buf);
        double w_xxx = d_xxx(ix,iy,iz,wp_buf);
        double p_xxx = d_xxx(ix,iy,iz,pp_buf);
        double u_yyy = d_yyy(ix,iy,iz,up_buf);
        double v_yyy = d_yyy(ix,iy,iz,vp_buf);
        double w_yyy = d_yyy(ix,iy,iz,wp_buf);
        double p_yyy = d_yyy(ix,iy,iz,pp_buf);
        double u_zzz = d_zzz(ix,iy,iz,up_buf);
        double v_zzz = d_zzz(ix,iy,iz,vp_buf);
        double w_zzz = d_zzz(ix,iy,iz,wp_buf);
        double p_zzz = d_zzz(ix,iy,iz,pp_buf);
        double u_xxy = d_xxy(ix,iy,iz,up_buf);
        double v_xxy = d_xxy(ix,iy,iz,vp_buf);
        double w_xxy = d_xxy(ix,iy,iz,wp_buf);
        double p_xxy = d_xxy(ix,iy,iz,pp_buf);
        double u_xxz = d_xxz(ix,iy,iz,up_buf);
        double v_xxz = d_xxz(ix,iy,iz,vp_buf);
        double w_xxz = d_xxz(ix,iy,iz,wp_buf);
        double p_xxz = d_xxz(ix,iy,iz,pp_buf);
        double u_xyy = d_xyy(ix,iy,iz,up_buf);
        double v_xyy = d_xyy(ix,iy,iz,vp_buf);
        double w_xyy = d_xyy(ix,iy,iz,wp_buf);
        double p_xyy = d_xyy(ix,iy,iz,pp_buf);
        double u_xzz = d_xzz(ix,iy,iz,up_buf);
        double v_xzz = d_xzz(ix,iy,iz,vp_buf);
        double w_xzz = d_xzz(ix,iy,iz,wp_buf);
        double p_xzz = d_xzz(ix,iy,iz,pp_buf);
        double u_yyz = d_yyz(ix,iy,iz,up_buf);
        double v_yyz = d_yyz(ix,iy,iz,vp_buf);
        double w_yyz = d_yyz(ix,iy,iz,wp_buf);
        double p_yyz = d_yyz(ix,iy,iz,pp_buf);
        double u_yzz = d_yzz(ix,iy,iz,up_buf);
        double v_yzz = d_yzz(ix,iy,iz,vp_buf);
        double w_yzz = d_yzz(ix,iy,iz,wp_buf);
        double p_yzz = d_yzz(ix,iy,iz,pp_buf);
        double u_xyz = d_xyz(ix,iy,iz,up_buf);
        double v_xyz = d_xyz(ix,iy,iz,vp_buf);
        double w_xyz = d_xyz(ix,iy,iz,wp_buf);

        double u_xxxx = d_xxxx(ix,iy,iz,up_buf);
        double u_xxxy = d_xxxy(ix,iy,iz,up_buf);
        double u_xxyy = d_xxyy(ix,iy,iz,up_buf);
        double u_xyyy = d_xyyy(ix,iy,iz,up_buf);
        double u_yyyy = d_yyyy(ix,iy,iz,up_buf);
        double u_xxxz = d_xxxz(ix,iy,iz,up_buf);
        double u_xyyz = d_xyyz(ix,iy,iz,up_buf);
        double u_xxzz = d_xxzz(ix,iy,iz,up_buf);
        double u_xyzz = d_xyzz(ix,iy,iz,up_buf);
        double u_yyzz = d_yyzz(ix,iy,iz,up_buf);
        double u_xzzz = d_xzzz(ix,iy,iz,up_buf);
        double u_zzzz = d_zzzz(ix,iy,iz,up_buf);

        double v_xxxx = d_xxxx(ix,iy,iz,vp_buf);
        double v_xxxy = d_xxxy(ix,iy,iz,vp_buf);
        double v_xxyy = d_xxyy(ix,iy,iz,vp_buf);
        double v_xyyy = d_xyyy(ix,iy,iz,vp_buf);
        double v_yyyy = d_yyyy(ix,iy,iz,vp_buf);
        double v_xxyz = d_xxyz(ix,iy,iz,vp_buf);
        double v_yyyz = d_yyyz(ix,iy,iz,vp_buf);
        double v_xxzz = d_xxzz(ix,iy,iz,vp_buf);
        double v_xyzz = d_xyzz(ix,iy,iz,vp_buf);
        double v_yyzz = d_yyzz(ix,iy,iz,vp_buf);
        double v_yzzz = d_yzzz(ix,iy,iz,vp_buf);
        double v_zzzz = d_zzzz(ix,iy,iz,vp_buf);

        double w_xxxx = d_xxxx(ix,iy,iz,wp_buf);
        double w_xxyy = d_xxyy(ix,iy,iz,wp_buf);
        double w_yyyy = d_yyyy(ix,iy,iz,wp_buf);
        double w_xxxz = d_xxxz(ix,iy,iz,wp_buf);
        double w_xxyz = d_xxyz(ix,iy,iz,wp_buf);
        double w_xyyz = d_xyyz(ix,iy,iz,wp_buf);
        double w_yyyz = d_yyyz(ix,iy,iz,wp_buf);
        double w_xxzz = d_xxzz(ix,iy,iz,wp_buf);
        double w_yyzz = d_yyzz(ix,iy,iz,wp_buf);
        double w_xzzz = d_xzzz(ix,iy,iz,wp_buf);
        double w_yzzz = d_yzzz(ix,iy,iz,wp_buf);
        double w_zzzz = d_zzzz(ix,iy,iz,wp_buf);

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

        rp_buf[ix-Ns][iy-Ns][iz-Ns] = rh_buf[ix][iy][iz] + 3*dt*r_t/2 + 5*dt*dt*r_tt/12;
        up_buf[ix-Ns][iy-Ns][iz-Ns] = uh_buf[ix][iy][iz] + 3*dt*u_t/2 + 5*dt*dt*u_tt/12;
        vp_buf[ix-Ns][iy-Ns][iz-Ns] = vh_buf[ix][iy][iz] + 3*dt*v_t/2 + 5*dt*dt*v_tt/12;
        wp_buf[ix-Ns][iy-Ns][iz-Ns] = wh_buf[ix][iy][iz] + 3*dt*w_t/2 + 5*dt*dt*w_tt/12;
        pp_buf[ix-Ns][iy-Ns][iz-Ns] = ph_buf[ix][iy][iz] + 3*dt*p_t/2 + 5*dt*dt*p_tt/12;

        rh_buf[ix-Ns][iy-Ns][iz-Ns] = rh_buf[ix][iy][iz] + dt*r_t;
        uh_buf[ix-Ns][iy-Ns][iz-Ns] = uh_buf[ix][iy][iz] + dt*u_t;
        vh_buf[ix-Ns][iy-Ns][iz-Ns] = vh_buf[ix][iy][iz] + dt*v_t;
        wh_buf[ix-Ns][iy-Ns][iz-Ns] = wh_buf[ix][iy][iz] + dt*w_t;
        ph_buf[ix-Ns][iy-Ns][iz-Ns] = ph_buf[ix][iy][iz] + dt*p_t;
      }
    }
  }
}

// Temporal BlockingでまとめてNTステップ更新
void update(param &p, state &s) {
  double rp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double up_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double vp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double wp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double pp_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];

  double rh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double uh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double vh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double wh_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];
  double ph_buf[NX+2*Ns][NY+2*Ns][NZ+2*Ns];

  // 床を読む
  int ix0 = p.x_origin;
  int iy0 = p.y_origin;
  int iz0 = p.z_origin;
  for(int ix = 0; ix < NX; ix++){
    for(int iy = 0; iy < NY; iy++){
      for(int iz = 0; iz < NZ; iz++){
        rp_buf[ix][iy][iz] = s.rp.get(ix0+ix, iy0+iy, iz0+iz);
        up_buf[ix][iy][iz] = s.up.get(ix0+ix, iy0+iy, iz0+iz);
        vp_buf[ix][iy][iz] = s.vp.get(ix0+ix, iy0+iy, iz0+iz);
        wp_buf[ix][iy][iz] = s.wp.get(ix0+ix, iy0+iy, iz0+iz);
        pp_buf[ix][iy][iz] = s.pp.get(ix0+ix, iy0+iy, iz0+iz);

        rh_buf[ix][iy][iz] = s.rh.get(ix0+ix, iy0+iy, iz0+iz);
        uh_buf[ix][iy][iz] = s.uh.get(ix0+ix, iy0+iy, iz0+iz);
        vh_buf[ix][iy][iz] = s.vh.get(ix0+ix, iy0+iy, iz0+iz);
        wh_buf[ix][iy][iz] = s.wh.get(ix0+ix, iy0+iy, iz0+iz);
        ph_buf[ix][iy][iz] = s.ph.get(ix0+ix, iy0+iy, iz0+iz);
      }
    }
  }

  for(int it = 0; it < NT; it++) {
    // x方向の壁を読む
    for(int ix = 0; ix < 2*Ns; ix++){
      for(int iy = 0; iy < NY+2*Ns; iy++){
        for(int iz = 0; iz < NZ+2*Ns; iz++){
          rp_buf[ix+NX][iy][iz] = p.rp_wall_x[it][ix][iy][iz];
          up_buf[ix+NX][iy][iz] = p.up_wall_x[it][ix][iy][iz];
          vp_buf[ix+NX][iy][iz] = p.vp_wall_x[it][ix][iy][iz];
          wp_buf[ix+NX][iy][iz] = p.wp_wall_x[it][ix][iy][iz];
          pp_buf[ix+NX][iy][iz] = p.pp_wall_x[it][ix][iy][iz];

          rh_buf[ix+NX][iy][iz] = p.rh_wall_x[it][ix][iy][iz];
          uh_buf[ix+NX][iy][iz] = p.uh_wall_x[it][ix][iy][iz];
          vh_buf[ix+NX][iy][iz] = p.vh_wall_x[it][ix][iy][iz];
          wh_buf[ix+NX][iy][iz] = p.wh_wall_x[it][ix][iy][iz];
          ph_buf[ix+NX][iy][iz] = p.ph_wall_x[it][ix][iy][iz];
        }
      }
    }
    // y方向の壁を読む
    for(int ix = 0; ix < NX+2*Ns; ix++){
      for(int iy = 0; iy < 2*Ns; iy++){
        for(int iz = 0; iz < NZ+2*Ns; iz++){
          rp_buf[ix][iy+NY][iz] = p.rp_wall_y[it][ix][iy][iz];
          up_buf[ix][iy+NY][iz] = p.up_wall_y[it][ix][iy][iz];
          vp_buf[ix][iy+NY][iz] = p.vp_wall_y[it][ix][iy][iz];
          wp_buf[ix][iy+NY][iz] = p.wp_wall_y[it][ix][iy][iz];
          pp_buf[ix][iy+NY][iz] = p.pp_wall_y[it][ix][iy][iz];

          rh_buf[ix][iy+NY][iz] = p.rh_wall_y[it][ix][iy][iz];
          uh_buf[ix][iy+NY][iz] = p.uh_wall_y[it][ix][iy][iz];
          vh_buf[ix][iy+NY][iz] = p.vh_wall_y[it][ix][iy][iz];
          wh_buf[ix][iy+NY][iz] = p.wh_wall_y[it][ix][iy][iz];
          ph_buf[ix][iy+NY][iz] = p.ph_wall_y[it][ix][iy][iz];
        }
      }
    }
    // z方向の壁を読む
    for(int ix = 0; ix < NX+2*Ns; ix++){
      for(int iy = 0; iy < NY+2*Ns; iy++){
        for(int iz = 0; iz < 2*Ns; iz++){
          rp_buf[ix][iy][iz+NZ] = p.rp_wall_z[it][ix][iy][iz];
          up_buf[ix][iy][iz+NZ] = p.up_wall_z[it][ix][iy][iz];
          vp_buf[ix][iy][iz+NZ] = p.vp_wall_z[it][ix][iy][iz];
          wp_buf[ix][iy][iz+NZ] = p.wp_wall_z[it][ix][iy][iz];
          pp_buf[ix][iy][iz+NZ] = p.pp_wall_z[it][ix][iy][iz];

          rh_buf[ix][iy][iz+NZ] = p.rh_wall_z[it][ix][iy][iz];
          uh_buf[ix][iy][iz+NZ] = p.uh_wall_z[it][ix][iy][iz];
          vh_buf[ix][iy][iz+NZ] = p.vh_wall_z[it][ix][iy][iz];
          wh_buf[ix][iy][iz+NZ] = p.wh_wall_z[it][ix][iy][iz];
          ph_buf[ix][iy][iz+NZ] = p.ph_wall_z[it][ix][iy][iz];
        }
      }
    }

    // x方向の壁を書く
    for(int ix = 0; ix < 2*Ns; ix++){
      for(int iy = 0; iy < NY+2*Ns; iy++){
        for(int iz = 0; iz < NZ+2*Ns; iz++){
          p.rp_wall_x[it][ix][iy][iz] = rp_buf[ix][iy][iz];
          p.up_wall_x[it][ix][iy][iz] = up_buf[ix][iy][iz];
          p.vp_wall_x[it][ix][iy][iz] = vp_buf[ix][iy][iz];
          p.wp_wall_x[it][ix][iy][iz] = wp_buf[ix][iy][iz];
          p.pp_wall_x[it][ix][iy][iz] = pp_buf[ix][iy][iz];

          p.rh_wall_x[it][ix][iy][iz] = rh_buf[ix][iy][iz];
          p.uh_wall_x[it][ix][iy][iz] = uh_buf[ix][iy][iz];
          p.vh_wall_x[it][ix][iy][iz] = vh_buf[ix][iy][iz];
          p.wh_wall_x[it][ix][iy][iz] = wh_buf[ix][iy][iz];
          p.ph_wall_x[it][ix][iy][iz] = ph_buf[ix][iy][iz];
        }
      }
    }
    // y方向の壁を書く
    for(int ix = 0; ix < NX+2*Ns; ix++){
      for(int iy = 0; iy < 2*Ns; iy++){
        for(int iz = 0; iz < NZ+2*Ns; iz++){
          p.rp_wall_y[it][ix][iy][iz] = rp_buf[ix][iy][iz];
          p.up_wall_y[it][ix][iy][iz] = up_buf[ix][iy][iz];
          p.vp_wall_y[it][ix][iy][iz] = vp_buf[ix][iy][iz];
          p.wp_wall_y[it][ix][iy][iz] = wp_buf[ix][iy][iz];
          p.pp_wall_y[it][ix][iy][iz] = pp_buf[ix][iy][iz];

          p.rh_wall_y[it][ix][iy][iz] = rh_buf[ix][iy][iz];
          p.uh_wall_y[it][ix][iy][iz] = uh_buf[ix][iy][iz];
          p.vh_wall_y[it][ix][iy][iz] = vh_buf[ix][iy][iz];
          p.wh_wall_y[it][ix][iy][iz] = wh_buf[ix][iy][iz];
          p.ph_wall_y[it][ix][iy][iz] = ph_buf[ix][iy][iz];
        }
      }
    }
    // z方向の壁を書く
    for(int ix = 0; ix < NX+2*Ns; ix++){
      for(int iy = 0; iy < NY+2*Ns; iy++){
        for(int iz = 0; iz < 2*Ns; iz++){
          p.rp_wall_z[it][ix][iy][iz] = rp_buf[ix][iy][iz];
          p.up_wall_z[it][ix][iy][iz] = up_buf[ix][iy][iz];
          p.vp_wall_z[it][ix][iy][iz] = vp_buf[ix][iy][iz];
          p.wp_wall_z[it][ix][iy][iz] = wp_buf[ix][iy][iz];
          p.pp_wall_z[it][ix][iy][iz] = pp_buf[ix][iy][iz];

          p.rh_wall_z[it][ix][iy][iz] = rh_buf[ix][iy][iz];
          p.uh_wall_z[it][ix][iy][iz] = uh_buf[ix][iy][iz];
          p.vh_wall_z[it][ix][iy][iz] = vh_buf[ix][iy][iz];
          p.wh_wall_z[it][ix][iy][iz] = wh_buf[ix][iy][iz];
          p.ph_wall_z[it][ix][iy][iz] = ph_buf[ix][iy][iz];
        }
      }
    }

    step(rp_buf,up_buf,vp_buf,wp_buf,pp_buf,rh_buf,uh_buf,vh_buf,wh_buf,ph_buf);
  }
  // 床を書く
  for(int ix = 0; ix < NX; ix++){
    for(int iy = 0; iy < NY; iy++){
      for(int iz = 0; iz < NZ; iz++){
        p.rp_res[ix][iy][iz] = rp_buf[ix][iy][iz];
        p.up_res[ix][iy][iz] = up_buf[ix][iy][iz];
        p.vp_res[ix][iy][iz] = vp_buf[ix][iy][iz];
        p.wp_res[ix][iy][iz] = wp_buf[ix][iy][iz];
        p.pp_res[ix][iy][iz] = pp_buf[ix][iy][iz];

        p.rh_res[ix][iy][iz] = rh_buf[ix][iy][iz];
        p.uh_res[ix][iy][iz] = uh_buf[ix][iy][iz];
        p.vh_res[ix][iy][iz] = vh_buf[ix][iy][iz];
        p.wh_res[ix][iy][iz] = wh_buf[ix][iy][iz];
        p.ph_res[ix][iy][iz] = ph_buf[ix][iy][iz];
      }
    }
  }
}

void fill_floor(state &s, int ix0, int ix1, int iy0, int iy1, int iz0, int iz1, int dx, int dy, int dz) {
  for(int ix = ix0; ix < ix1; ix++){
    for(int iy = iy0; iy < iy1; iy++){
      for(int iz = iz0; iz < iz1; iz++){
        s.rp.get(ix, iy, iz) = s.rp.get(ix+dx, iy+dy, iz+dz);
        s.up.get(ix, iy, iz) = s.up.get(ix+dx, iy+dy, iz+dz);
        s.vp.get(ix, iy, iz) = s.vp.get(ix+dx, iy+dy, iz+dz);
        s.wp.get(ix, iy, iz) = s.wp.get(ix+dx, iy+dy, iz+dz);
        s.pp.get(ix, iy, iz) = s.pp.get(ix+dx, iy+dy, iz+dz);

        s.rh.get(ix, iy, iz) = s.rh.get(ix+dx, iy+dy, iz+dz);
        s.uh.get(ix, iy, iz) = s.uh.get(ix+dx, iy+dy, iz+dz);
        s.vh.get(ix, iy, iz) = s.vh.get(ix+dx, iy+dy, iz+dz);
        s.wh.get(ix, iy, iz) = s.wh.get(ix+dx, iy+dy, iz+dz);
        s.ph.get(ix, iy, iz) = s.ph.get(ix+dx, iy+dy, iz+dz);
      }
    }
  }
}

// 壁出力を一時的に保持する配列
double rp_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double up_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double vp_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double wp_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double pp_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double rh_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double uh_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double vh_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double wh_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];
double ph_tmp_wall_x[MY+2][MZ+2][NT][2*Ns][NY+2*Ns][NZ+2*Ns];

double rp_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double up_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double vp_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double wp_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double pp_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double rh_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double uh_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double vh_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double wh_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];
double ph_tmp_wall_y[MZ+2][NT][NX+2*Ns][2*Ns][NZ+2*Ns];

struct tmp_buf {
  // 領域を分割し、updateを呼び全体を時間発展させる
  // 床出力を一時的に保持する配列
  double rp_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double up_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double vp_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double wp_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double pp_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double rh_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double uh_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double vh_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double wh_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
  double ph_tmp[(MX+2)*NX][(MY+2)*NY][(MZ+2)*NZ];
};

void next(navi &n, state &s) {
  auto tmp = unique_ptr<tmp_buf>(new tmp_buf);

  param p;

  // 床の準備(本来はMPI通信が必要)
  fill_floor(s, 0,2*NX,0,2*NY,0,2*NZ,MX*NX,MY*NY,MZ*NZ);
  fill_floor(s, 2*NX,(MX+2)*NX,0,2*NY,0,2*NZ,0,MY*NY,MZ*NZ);
  fill_floor(s, 0,2*NX,2*NY,(MY+2)*NY,0,2*NZ,MX*NX,0,MZ*NZ);
  fill_floor(s, 0,2*NX,0,2*NY,2*NZ,(MZ+2)*NZ,MX*NX,MY*NY,0);
  fill_floor(s, 0,2*NX,2*NY,(MY+2)*NY,2*NZ,(MZ+2)*NZ,MX*NX,0,0);
  fill_floor(s, 2*NX,(MX+2)*NX,0,2*NY,2*NZ,(MZ+2)*NZ,0,MY*NY,0);
  fill_floor(s, 2*NX,(MX+2)*NX,2*NY,(MY+2)*NY,0,2*NZ,0,0,MZ*NZ);

  // tmp_wall_xの初期化
  for (int jz = 0; jz < MZ+2; jz++) {
    for (int jy = 0; jy < MY+2; jy++) {
      for(int it = 0; it < NT; it++) {
        for(int ix = 0; ix < 2*Ns; ix++){
          for(int iy = 0; iy < NY+2*Ns; iy++){
            for(int iz = 0; iz < NZ+2*Ns; iz++){
              rp_tmp_wall_x[jy][jz][it][ix][iy][iz] = 1.0;
              up_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              vp_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              wp_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              pp_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;

              rh_tmp_wall_x[jy][jz][it][ix][iy][iz] = 1.0;
              uh_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              vh_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              wh_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
              ph_tmp_wall_x[jy][jz][it][ix][iy][iz] = 0.0;
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
              rp_tmp_wall_y[jz][it][ix][iy][iz] = 1.0;
              up_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              vp_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              wp_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              pp_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;

              rh_tmp_wall_y[jz][it][ix][iy][iz] = 1.0;
              uh_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              vh_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              wh_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
              ph_tmp_wall_y[jz][it][ix][iy][iz] = 0.0;
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
              p.rp_wall_z[it][ix][iy][iz] = 1.0;
              p.up_wall_z[it][ix][iy][iz] = 0.0;
              p.vp_wall_z[it][ix][iy][iz] = 0.0;
              p.wp_wall_z[it][ix][iy][iz] = 0.0;
              p.pp_wall_z[it][ix][iy][iz] = 0.0;

              p.rh_wall_z[it][ix][iy][iz] = 1.0;
              p.uh_wall_z[it][ix][iy][iz] = 0.0;
              p.vh_wall_z[it][ix][iy][iz] = 0.0;
              p.wh_wall_z[it][ix][iy][iz] = 0.0;
              p.ph_wall_z[it][ix][iy][iz] = 0.0;
            }
          }
        }
      }

      for(int jz = MZ+1;jz >= 0; jz--) {
        p.x_origin = jx*NX;
        p.y_origin = jy*NY;
        p.z_origin = jz*NZ;

        for(int it = 0; it < NT; it++) {
          // x壁を渡す
          for(int ix = 0; ix < 2*Ns; ix++){
            for(int iy = 0; iy < NY+2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p.rp_wall_x[it][ix][iy][iz] = rp_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.up_wall_x[it][ix][iy][iz] = up_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.vp_wall_x[it][ix][iy][iz] = vp_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.wp_wall_x[it][ix][iy][iz] = wp_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.pp_wall_x[it][ix][iy][iz] = pp_tmp_wall_x[jy][jz][it][ix][iy][iz];

                p.rh_wall_x[it][ix][iy][iz] = rh_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.uh_wall_x[it][ix][iy][iz] = uh_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.vh_wall_x[it][ix][iy][iz] = vh_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.wh_wall_x[it][ix][iy][iz] = wh_tmp_wall_x[jy][jz][it][ix][iy][iz];
                p.ph_wall_x[it][ix][iy][iz] = ph_tmp_wall_x[jy][jz][it][ix][iy][iz];
              }
            }
          }

          // y壁を渡す
          for(int ix = 0; ix < NX+2*Ns; ix++){
            for(int iy = 0; iy < 2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                p.rp_wall_y[it][ix][iy][iz] = rp_tmp_wall_y[jz][it][ix][iy][iz];
                p.up_wall_y[it][ix][iy][iz] = up_tmp_wall_y[jz][it][ix][iy][iz];
                p.vp_wall_y[it][ix][iy][iz] = vp_tmp_wall_y[jz][it][ix][iy][iz];
                p.wp_wall_y[it][ix][iy][iz] = wp_tmp_wall_y[jz][it][ix][iy][iz];
                p.pp_wall_y[it][ix][iy][iz] = pp_tmp_wall_y[jz][it][ix][iy][iz];

                p.rh_wall_y[it][ix][iy][iz] = rh_tmp_wall_y[jz][it][ix][iy][iz];
                p.uh_wall_y[it][ix][iy][iz] = uh_tmp_wall_y[jz][it][ix][iy][iz];
                p.vh_wall_y[it][ix][iy][iz] = vh_tmp_wall_y[jz][it][ix][iy][iz];
                p.wh_wall_y[it][ix][iy][iz] = wh_tmp_wall_y[jz][it][ix][iy][iz];
                p.ph_wall_y[it][ix][iy][iz] = ph_tmp_wall_y[jz][it][ix][iy][iz];
              }
            }
          }
        }

        update(p, s);

        for(int it = 0; it < NT; it++) {
          // x壁の受けとり
          for(int ix = 0; ix < 2*Ns; ix++){
            for(int iy = 0; iy < NY+2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                rp_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.rp_wall_x[it][ix][iy][iz];
                up_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.up_wall_x[it][ix][iy][iz];
                vp_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.vp_wall_x[it][ix][iy][iz];
                wp_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.wp_wall_x[it][ix][iy][iz];
                pp_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.pp_wall_x[it][ix][iy][iz];

                rh_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.rh_wall_x[it][ix][iy][iz];
                uh_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.uh_wall_x[it][ix][iy][iz];
                vh_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.vh_wall_x[it][ix][iy][iz];
                wh_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.wh_wall_x[it][ix][iy][iz];
                ph_tmp_wall_x[jy][jz][it][ix][iy][iz] = p.ph_wall_x[it][ix][iy][iz];
              }
            }
          }

          // y壁の受けとり
          for(int ix = 0; ix < NX+2*Ns; ix++){
            for(int iy = 0; iy < 2*Ns; iy++){
              for(int iz = 0; iz < NZ+2*Ns; iz++){
                rp_tmp_wall_y[jz][it][ix][iy][iz] = p.rp_wall_y[it][ix][iy][iz];
                up_tmp_wall_y[jz][it][ix][iy][iz] = p.up_wall_y[it][ix][iy][iz];
                vp_tmp_wall_y[jz][it][ix][iy][iz] = p.vp_wall_y[it][ix][iy][iz];
                wp_tmp_wall_y[jz][it][ix][iy][iz] = p.wp_wall_y[it][ix][iy][iz];
                pp_tmp_wall_y[jz][it][ix][iy][iz] = p.pp_wall_y[it][ix][iy][iz];

                rh_tmp_wall_y[jz][it][ix][iy][iz] = p.rh_wall_y[it][ix][iy][iz];
                uh_tmp_wall_y[jz][it][ix][iy][iz] = p.uh_wall_y[it][ix][iy][iz];
                vh_tmp_wall_y[jz][it][ix][iy][iz] = p.vh_wall_y[it][ix][iy][iz];
                wh_tmp_wall_y[jz][it][ix][iy][iz] = p.wh_wall_y[it][ix][iy][iz];
                ph_tmp_wall_y[jz][it][ix][iy][iz] = p.ph_wall_y[it][ix][iy][iz];
              }
            }
          }
        }

        // 床の受けとり
        for(int ix = 0; ix < NX; ix++) {
          for(int iy = 0; iy < NY; iy++) {
            for(int iz = 0; iz < NZ; iz++) {
              tmp->rp_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.rp_res[ix][iy][iz];
              tmp->up_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.up_res[ix][iy][iz];
              tmp->vp_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.vp_res[ix][iy][iz];
              tmp->wp_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.wp_res[ix][iy][iz];
              tmp->pp_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.pp_res[ix][iy][iz];

              tmp->rh_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.rh_res[ix][iy][iz];
              tmp->uh_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.uh_res[ix][iy][iz];
              tmp->vh_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.vh_res[ix][iy][iz];
              tmp->wh_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.wh_res[ix][iy][iz];
              tmp->ph_tmp[ix+jx*NX][iy+jy*NY][iz+jz*NZ] = p.ph_res[ix][iy][iz];
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
        s.rp.get(x, y, z) = tmp->rp_tmp[ix][iy][iz];
        s.up.get(x, y, z) = tmp->up_tmp[ix][iy][iz];
        s.vp.get(x, y, z) = tmp->vp_tmp[ix][iy][iz];
        s.wp.get(x, y, z) = tmp->wp_tmp[ix][iy][iz];
        s.pp.get(x, y, z) = tmp->pp_tmp[ix][iy][iz];

        s.rh.get(x, y, z) = tmp->rh_tmp[ix][iy][iz];
        s.uh.get(x, y, z) = tmp->uh_tmp[ix][iy][iz];
        s.vh.get(x, y, z) = tmp->vh_tmp[ix][iy][iz];
        s.wh.get(x, y, z) = tmp->wh_tmp[ix][iy][iz];
        s.ph.get(x, y, z) = tmp->ph_tmp[ix][iy][iz];
      }
    }
  }

  // NTだけタイムステップを進める
  n.time_step += NT;
}
