vis1 = 4.0*pow(3.0,-1.0)*u_xx + pow(3.0,-1.0)*v_xy + pow(3.0,-1.0)*w_xz + u_yy + u_zz
vis2 = pow(3.0,-1.0)*u_xy + 4.0*pow(3.0,-1.0)*v_yy + pow(3.0,-1.0)*w_yz + v_xx + v_zz
vis3 = pow(3.0,-1.0)*u_xz + pow(3.0,-1.0)*v_yz + 4.0*pow(3.0,-1.0)*w_zz + w_xx + w_yy
vis1_x = 4.0*pow(3.0,-1.0)*u_xxx + pow(3.0,-1.0)*v_xxy + pow(3.0,-1.0)*w_xxz + u_xyy + u_xzz
vis2_x = pow(3.0,-1.0)*u_xxy + 4.0*pow(3.0,-1.0)*v_xyy + pow(3.0,-1.0)*w_xyz + v_xzz + v_xxx
vis3_x = pow(3.0,-1.0)*u_xxz + pow(3.0,-1.0)*v_xyz + 4.0*pow(3.0,-1.0)*w_xzz + w_xyy + w_xxx
vis1_y = 4.0*pow(3.0,-1.0)*u_xxy + pow(3.0,-1.0)*v_xyy + pow(3.0,-1.0)*w_xyz + u_yzz + u_yyy
vis2_y = pow(3.0,-1.0)*u_xyy + 4.0*pow(3.0,-1.0)*v_yyy + pow(3.0,-1.0)*w_yyz + v_xxy + v_yzz
vis3_y = pow(3.0,-1.0)*u_xyz + pow(3.0,-1.0)*v_yyz + 4.0*pow(3.0,-1.0)*w_yzz + w_xxy + w_yyy
vis1_z = 4.0*pow(3.0,-1.0)*u_xxz + pow(3.0,-1.0)*v_xyz + pow(3.0,-1.0)*w_xzz + u_yyz + u_zzz
vis2_z = pow(3.0,-1.0)*u_xyz + 4.0*pow(3.0,-1.0)*v_yyz + pow(3.0,-1.0)*w_yzz + v_xxz + v_zzz
vis3_z = pow(3.0,-1.0)*u_xzz + pow(3.0,-1.0)*v_yzz + 4.0*pow(3.0,-1.0)*w_zzz + w_xxz + w_yyz
vis1_xx = 4.0*pow(3.0,-1.0)*u_xxxx + pow(3.0,-1.0)*v_xxxy + pow(3.0,-1.0)*w_xxxz + u_xxyy + u_xxzz
vis2_xx = pow(3.0,-1.0)*u_xxxy + 4.0*pow(3.0,-1.0)*v_xxyy + pow(3.0,-1.0)*w_xxyz + v_xxzz + v_xxxx
vis3_xx = pow(3.0,-1.0)*u_xxxz + pow(3.0,-1.0)*v_xxyz + 4.0*pow(3.0,-1.0)*w_xxzz + w_xxyy + w_xxxx
vis1_yy = 4.0*pow(3.0,-1.0)*u_xxyy + pow(3.0,-1.0)*v_xyyy + pow(3.0,-1.0)*w_xyyz + u_yyzz + u_yyyy
vis2_yy = pow(3.0,-1.0)*u_xyyy + 4.0*pow(3.0,-1.0)*v_yyyy + pow(3.0,-1.0)*w_yyyz + v_xxyy + v_yyzz
vis3_yy = pow(3.0,-1.0)*u_xyyz + pow(3.0,-1.0)*v_yyyz + 4.0*pow(3.0,-1.0)*w_yyzz + w_xxyy + w_yyyy
vis1_zz = 4.0*pow(3.0,-1.0)*u_xxzz + pow(3.0,-1.0)*v_xyzz + pow(3.0,-1.0)*w_xzzz + u_yyzz + u_zzzz
vis2_zz = pow(3.0,-1.0)*u_xyzz + 4.0*pow(3.0,-1.0)*v_yyzz + pow(3.0,-1.0)*w_yzzz + v_xxzz + v_zzzz
vis3_zz = pow(3.0,-1.0)*u_xzzz + pow(3.0,-1.0)*v_yzzz + 4.0*pow(3.0,-1.0)*w_zzzz + w_xxzz + w_yyzz
vis1_xy = 4.0*pow(3.0,-1.0)*u_xxxy + pow(3.0,-1.0)*v_xxyy + pow(3.0,-1.0)*w_xxyz + u_xyzz + u_xyyy
vis2_xy = pow(3.0,-1.0)*u_xxyy + 4.0*pow(3.0,-1.0)*v_xyyy + pow(3.0,-1.0)*w_xyyz + v_xyzz + v_xxxy
vis3_xy = pow(3.0,-1.0)*u_xxyz + pow(3.0,-1.0)*v_xyyz + 4.0*pow(3.0,-1.0)*w_xyzz + w_xyyy + w_xxxy
vis1_yz = 4.0*pow(3.0,-1.0)*u_xxyz + pow(3.0,-1.0)*v_xyyz + pow(3.0,-1.0)*w_xyzz + u_yzzz + u_yyyz
vis2_yz = pow(3.0,-1.0)*u_xyyz + 4.0*pow(3.0,-1.0)*v_yyyz + pow(3.0,-1.0)*w_yyzz + v_xxyz + v_yzzz
vis3_yz = pow(3.0,-1.0)*u_xyzz + pow(3.0,-1.0)*v_yyzz + 4.0*pow(3.0,-1.0)*w_yzzz + w_xxyz + w_yyyz
vis1_xz = 4.0*pow(3.0,-1.0)*u_xxxz + pow(3.0,-1.0)*v_xxyz + pow(3.0,-1.0)*w_xxzz + u_xyyz + u_xzzz
vis2_xz = pow(3.0,-1.0)*u_xxyz + 4.0*pow(3.0,-1.0)*v_xyyz + pow(3.0,-1.0)*w_xyzz + v_xzzz + v_xxxz
vis3_xz = pow(3.0,-1.0)*u_xxzz + pow(3.0,-1.0)*v_xyzz + 4.0*pow(3.0,-1.0)*w_xzzz + w_xyyz + w_xxxz
vis1_t = 4.0*pow(3.0,-1.0)*u_txx + pow(3.0,-1.0)*v_txy + pow(3.0,-1.0)*w_txz + u_tyy + u_tzz
vis2_t = pow(3.0,-1.0)*u_txy + 4.0*pow(3.0,-1.0)*v_tyy + pow(3.0,-1.0)*w_tyz + v_txx + v_tzz
vis3_t = pow(3.0,-1.0)*u_txz + pow(3.0,-1.0)*v_tyz + 4.0*pow(3.0,-1.0)*w_tzz + w_txx + w_tyy

r_t = -r*u_x - r*v_y - r*w_z - r_x*u - r_y*v - r_z*w
u_t = c*pow(r,-1.0)*vis1 - p_x*pow(r,-1.0) - u*u_x - u_y*v - u_z*w
v_t = c*pow(r,-1.0)*vis2 - p_y*pow(r,-1.0) - u*v_x - v*v_y - v_z*w
w_t = c*pow(r,-1.0)*vis3 - p_z*pow(r,-1.0) - u*w_x - v*w_y - w*w_z
p_t = -c2*u*vis1 - c2*v*vis2 - c2*vis3*w - gm*p*u_x - gm*p*v_y - gm*p*w_z - p_x*u - p_y*v - p_z*w

r_tx = -r*u_xx - r*v_xy - r*w_xz - 2.0*r_x*u_x - r_x*v_y - r_x*w_z - r_xy*v - r_xz*w - r_xx*u - r_y*v_x - r_z*w_x
u_tx = -c*pow(r,-2.0)*r_x*vis1 + c*pow(r,-1.0)*vis1_x + p_x*pow(r,-2.0)*r_x - p_xx*pow(r,-1.0) - u*u_xx - pow(u_x,2.0) - u_xy*v - u_xz*w - u_y*v_x - u_z*w_x
v_tx = -c*pow(r,-2.0)*r_x*vis2 + c*pow(r,-1.0)*vis2_x - p_xy*pow(r,-1.0) + p_y*pow(r,-2.0)*r_x - u*v_xx - u_x*v_x - v*v_xy - v_x*v_y - v_xz*w - v_z*w_x
w_tx = -c*pow(r,-2.0)*r_x*vis3 + c*pow(r,-1.0)*vis3_x - p_xz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_x - u*w_xx - u_x*w_x - v*w_xy - v_x*w_y - w*w_xz - w_x*w_z
p_tx = -c2*u*vis1_x - c2*u_x*vis1 - c2*v*vis2_x - c2*v_x*vis2 - c2*vis3*w_x - c2*vis3_x*w - gm*p*u_xx - gm*p*v_xy - gm*p*w_xz - gm*p_x*u_x - gm*p_x*v_y - gm*p_x*w_z - p_x*u_x - p_xy*v - p_xz*w - p_xx*u - p_y*v_x - p_z*w_x

r_ty = -r*u_xy - r*v_yy - r*w_yz - r_x*u_y - r_xy*u - r_y*u_x - 2.0*r_y*v_y - r_y*w_z - r_yz*w - r_yy*v - r_z*w_y
u_ty = -c*pow(r,-2.0)*r_y*vis1 + c*pow(r,-1.0)*vis1_y + p_x*pow(r,-2.0)*r_y - p_xy*pow(r,-1.0) - u*u_xy - u_x*u_y - u_y*v_y - u_yz*w - u_yy*v - u_z*w_y
v_ty = -c*pow(r,-2.0)*r_y*vis2 + c*pow(r,-1.0)*vis2_y + p_y*pow(r,-2.0)*r_y - p_yy*pow(r,-1.0) - u*v_xy - u_y*v_x - v*v_yy - pow(v_y,2.0) - v_yz*w - v_z*w_y
w_ty = -c*pow(r,-2.0)*r_y*vis3 + c*pow(r,-1.0)*vis3_y - p_yz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_y - u*w_xy - u_y*w_x - v*w_yy - v_y*w_y - w*w_yz - w_y*w_z
p_ty = -c2*u*vis1_y - c2*u_y*vis1 - c2*v*vis2_y - c2*v_y*vis2 - c2*vis3*w_y - c2*vis3_y*w - gm*p*u_xy - gm*p*v_yy - gm*p*w_yz - gm*p_y*u_x - gm*p_y*v_y - gm*p_y*w_z - p_x*u_y - p_xy*u - p_y*v_y - p_yz*w - p_yy*v - p_z*w_y

r_tz = -r*u_xz - r*v_yz - r*w_zz - r_x*u_z - r_xz*u - r_y*v_z - r_yz*v - r_z*u_x - r_z*v_y - 2.0*r_z*w_z - r_zz*w
u_tz = -c*pow(r,-2.0)*r_z*vis1 + c*pow(r,-1.0)*vis1_z + p_x*pow(r,-2.0)*r_z - p_xz*pow(r,-1.0) - u*u_xz - u_x*u_z - u_y*v_z - u_yz*v - u_z*w_z - u_zz*w
v_tz = -c*pow(r,-2.0)*r_z*vis2 + c*pow(r,-1.0)*vis2_z + p_y*pow(r,-2.0)*r_z - p_yz*pow(r,-1.0) - u*v_xz - u_z*v_x - v*v_yz - v_y*v_z - v_z*w_z - v_zz*w
w_tz = -c*pow(r,-2.0)*r_z*vis3 + c*pow(r,-1.0)*vis3_z + p_z*pow(r,-2.0)*r_z - p_zz*pow(r,-1.0) - u*w_xz - u_z*w_x - v*w_yz - v_z*w_y - w*w_zz - pow(w_z,2.0)
p_tz = -c2*u*vis1_z - c2*u_z*vis1 - c2*v*vis2_z - c2*v_z*vis2 - c2*vis3*w_z - c2*vis3_z*w - gm*p*u_xz - gm*p*v_yz - gm*p*w_zz - gm*p_z*u_x - gm*p_z*v_y - gm*p_z*w_z - p_x*u_z - p_xz*u - p_y*v_z - p_yz*v - p_z*w_z - p_zz*w

r_txx = -r*u_xxx - r*v_xxy - r*w_xxz - 3.0*r_x*u_xx - 2.0*r_x*v_xy - 2.0*r_x*w_xz - 2.0*r_xy*v_x - 2.0*r_xz*w_x - 3.0*r_xx*u_x - r_xx*v_y - r_xx*w_z - r_xxy*v - r_xxz*w - r_xxx*u - r_y*v_xx - r_z*w_xx
u_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_x*vis1_x - c*pow(r,-2.0)*r_xx*vis1 + c*pow(r,-1.0)*vis1_xx - 2.0*p_x*pow(r,-3.0)*pow(r_x,2.0) + p_x*pow(r,-2.0)*r_xx + 2.0*p_xx*pow(r,-2.0)*r_x - p_xxx*pow(r,-1.0) - u*u_xxx - 3.0*u_x*u_xx - 2.0*u_xy*v_x - 2.0*u_xz*w_x - u_xxy*v - u_xxz*w - u_y*v_xx - u_z*w_xx
v_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_x*vis2_x - c*pow(r,-2.0)*r_xx*vis2 + c*pow(r,-1.0)*vis2_xx + 2.0*p_xy*pow(r,-2.0)*r_x - p_xxy*pow(r,-1.0) - 2.0*p_y*pow(r,-3.0)*pow(r_x,2.0) + p_y*pow(r,-2.0)*r_xx - u*v_xxx - 2.0*u_x*v_xx - u_xx*v_x - v*v_xxy - 2.0*v_x*v_xy - 2.0*v_xz*w_x - v_xx*v_y - v_xxz*w - v_z*w_xx
w_txx = 2.0*c*pow(r,-3.0)*pow(r_x,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_x*vis3_x - c*pow(r,-2.0)*r_xx*vis3 + c*pow(r,-1.0)*vis3_xx + 2.0*p_xz*pow(r,-2.0)*r_x - p_xxz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*pow(r_x,2.0) + p_z*pow(r,-2.0)*r_xx - u*w_xxx - 2.0*u_x*w_xx - u_xx*w_x - v*w_xxy - 2.0*v_x*w_xy - v_xx*w_y - w*w_xxz - 2.0*w_x*w_xz - w_xx*w_z
p_txx = -c2*u*vis1_xx - 2.0*c2*u_x*vis1_x - c2*u_xx*vis1 - c2*v*vis2_xx - 2.0*c2*v_x*vis2_x - c2*v_xx*vis2 - c2*vis3*w_xx - 2.0*c2*vis3_x*w_x - c2*vis3_xx*w - gm*p*u_xxx - gm*p*v_xxy - gm*p*w_xxz - 2.0*gm*p_x*u_xx - 2.0*gm*p_x*v_xy - 2.0*gm*p_x*w_xz - gm*p_xx*u_x - gm*p_xx*v_y - gm*p_xx*w_z - p_x*u_xx - 2.0*p_xy*v_x - 2.0*p_xz*w_x - 2.0*p_xx*u_x - p_xxy*v - p_xxz*w - p_xxx*u - p_y*v_xx - p_z*w_xx

r_tyy = -r*u_xyy - r*v_yyy - r*w_yyz - r_x*u_yy - 2.0*r_xy*u_y - r_xyy*u - 2.0*r_y*u_xy - 3.0*r_y*v_yy - 2.0*r_y*w_yz - 2.0*r_yz*w_y - r_yy*u_x - 3.0*r_yy*v_y - r_yy*w_z - r_yyz*w - r_yyy*v - r_z*w_yy
u_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_y*vis1_y - c*pow(r,-2.0)*r_yy*vis1 + c*pow(r,-1.0)*vis1_yy - 2.0*p_x*pow(r,-3.0)*pow(r_y,2.0) + p_x*pow(r,-2.0)*r_yy + 2.0*p_xy*pow(r,-2.0)*r_y - p_xyy*pow(r,-1.0) - u*u_xyy - u_x*u_yy - 2.0*u_xy*u_y - u_y*v_yy - 2.0*u_yz*w_y - 2.0*u_yy*v_y - u_yyz*w - u_yyy*v - u_z*w_yy
v_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_y*vis2_y - c*pow(r,-2.0)*r_yy*vis2 + c*pow(r,-1.0)*vis2_yy - 2.0*p_y*pow(r,-3.0)*pow(r_y,2.0) + p_y*pow(r,-2.0)*r_yy + 2.0*p_yy*pow(r,-2.0)*r_y - p_yyy*pow(r,-1.0) - u*v_xyy - 2.0*u_y*v_xy - u_yy*v_x - v*v_yyy - 3.0*v_y*v_yy - 2.0*v_yz*w_y - v_yyz*w - v_z*w_yy
w_tyy = 2.0*c*pow(r,-3.0)*pow(r_y,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_y*vis3_y - c*pow(r,-2.0)*r_yy*vis3 + c*pow(r,-1.0)*vis3_yy + 2.0*p_yz*pow(r,-2.0)*r_y - p_yyz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*pow(r_y,2.0) + p_z*pow(r,-2.0)*r_yy - u*w_xyy - 2.0*u_y*w_xy - u_yy*w_x - v*w_yyy - 2.0*v_y*w_yy - v_yy*w_y - w*w_yyz - 2.0*w_y*w_yz - w_yy*w_z
p_tyy = -c2*u*vis1_yy - 2.0*c2*u_y*vis1_y - c2*u_yy*vis1 - c2*v*vis2_yy - 2.0*c2*v_y*vis2_y - c2*v_yy*vis2 - c2*vis3*w_yy - 2.0*c2*vis3_y*w_y - c2*vis3_yy*w - gm*p*u_xyy - gm*p*v_yyy - gm*p*w_yyz - 2.0*gm*p_y*u_xy - 2.0*gm*p_y*v_yy - 2.0*gm*p_y*w_yz - gm*p_yy*u_x - gm*p_yy*v_y - gm*p_yy*w_z - p_x*u_yy - 2.0*p_xy*u_y - p_xyy*u - p_y*v_yy - 2.0*p_yz*w_y - 2.0*p_yy*v_y - p_yyz*w - p_yyy*v - p_z*w_yy

r_tzz = -r*u_xzz - r*v_yzz - r*w_zzz - r_x*u_zz - 2.0*r_xz*u_z - r_xzz*u - r_y*v_zz - 2.0*r_yz*v_z - r_yzz*v - 2.0*r_z*u_xz - 2.0*r_z*v_yz - 3.0*r_z*w_zz - r_zz*u_x - r_zz*v_y - 3.0*r_zz*w_z - r_zzz*w
u_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis1 - 2.0*c*pow(r,-2.0)*r_z*vis1_z - c*pow(r,-2.0)*r_zz*vis1 + c*pow(r,-1.0)*vis1_zz - 2.0*p_x*pow(r,-3.0)*pow(r_z,2.0) + p_x*pow(r,-2.0)*r_zz + 2.0*p_xz*pow(r,-2.0)*r_z - p_xzz*pow(r,-1.0) - u*u_xzz - u_x*u_zz - 2.0*u_xz*u_z - u_y*v_zz - 2.0*u_yz*v_z - u_yzz*v - u_z*w_zz - 2.0*u_zz*w_z - u_zzz*w
v_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis2 - 2.0*c*pow(r,-2.0)*r_z*vis2_z - c*pow(r,-2.0)*r_zz*vis2 + c*pow(r,-1.0)*vis2_zz - 2.0*p_y*pow(r,-3.0)*pow(r_z,2.0) + p_y*pow(r,-2.0)*r_zz + 2.0*p_yz*pow(r,-2.0)*r_z - p_yzz*pow(r,-1.0) - u*v_xzz - 2.0*u_z*v_xz - u_zz*v_x - v*v_yzz - v_y*v_zz - 2.0*v_yz*v_z - v_z*w_zz - 2.0*v_zz*w_z - v_zzz*w
w_tzz = 2.0*c*pow(r,-3.0)*pow(r_z,2.0)*vis3 - 2.0*c*pow(r,-2.0)*r_z*vis3_z - c*pow(r,-2.0)*r_zz*vis3 + c*pow(r,-1.0)*vis3_zz - 2.0*p_z*pow(r,-3.0)*pow(r_z,2.0) + p_z*pow(r,-2.0)*r_zz + 2.0*p_zz*pow(r,-2.0)*r_z - p_zzz*pow(r,-1.0) - u*w_xzz - 2.0*u_z*w_xz - u_zz*w_x - v*w_yzz - 2.0*v_z*w_yz - v_zz*w_y - w*w_zzz - 3.0*w_z*w_zz
p_tzz = -c2*u*vis1_zz - 2.0*c2*u_z*vis1_z - c2*u_zz*vis1 - c2*v*vis2_zz - 2.0*c2*v_z*vis2_z - c2*v_zz*vis2 - c2*vis3*w_zz - 2.0*c2*vis3_z*w_z - c2*vis3_zz*w - gm*p*u_xzz - gm*p*v_yzz - gm*p*w_zzz - 2.0*gm*p_z*u_xz - 2.0*gm*p_z*v_yz - 2.0*gm*p_z*w_zz - gm*p_zz*u_x - gm*p_zz*v_y - gm*p_zz*w_z - p_x*u_zz - 2.0*p_xz*u_z - p_xzz*u - p_y*v_zz - 2.0*p_yz*v_z - p_yzz*v - p_z*w_zz - 2.0*p_zz*w_z - p_zzz*w

r_txy = -r*u_xxy - r*v_xyy - r*w_xyz - 2.0*r_x*u_xy - r_x*v_yy - r_x*w_yz - 2.0*r_xy*u_x - 2.0*r_xy*v_y - r_xy*w_z - r_xyz*w - r_xyy*v - r_xz*w_y - r_xx*u_y - r_xxy*u - r_y*u_xx - 2.0*r_y*v_xy - r_y*w_xz - r_yz*w_x - r_yy*v_x - r_z*w_xy
u_txy = 2.0*c*pow(r,-3.0)*r_x*r_y*vis1 - c*pow(r,-2.0)*r_x*vis1_y - c*pow(r,-2.0)*r_xy*vis1 - c*pow(r,-2.0)*r_y*vis1_x + c*pow(r,-1.0)*vis1_xy - 2.0*p_x*pow(r,-3.0)*r_x*r_y + p_x*pow(r,-2.0)*r_xy + p_xy*pow(r,-2.0)*r_x + p_xx*pow(r,-2.0)*r_y - p_xxy*pow(r,-1.0) - u*u_xxy - 2.0*u_x*u_xy - u_xy*v_y - u_xyz*w - u_xyy*v - u_xz*w_y - u_xx*u_y - u_y*v_xy - u_yz*w_x - u_yy*v_x - u_z*w_xy
v_txy = 2.0*c*pow(r,-3.0)*r_x*r_y*vis2 - c*pow(r,-2.0)*r_x*vis2_y - c*pow(r,-2.0)*r_xy*vis2 - c*pow(r,-2.0)*r_y*vis2_x + c*pow(r,-1.0)*vis2_xy + p_xy*pow(r,-2.0)*r_y - p_xyy*pow(r,-1.0) - 2.0*p_y*pow(r,-3.0)*r_x*r_y + p_y*pow(r,-2.0)*r_xy + p_yy*pow(r,-2.0)*r_x - u*v_xxy - u_x*v_xy - u_xy*v_x - u_y*v_xx - v*v_xyy - v_x*v_yy - 2.0*v_xy*v_y - v_xyz*w - v_xz*w_y - v_yz*w_x - v_z*w_xy
w_txy = 2.0*c*pow(r,-3.0)*r_x*r_y*vis3 - c*pow(r,-2.0)*r_x*vis3_y - c*pow(r,-2.0)*r_xy*vis3 - c*pow(r,-2.0)*r_y*vis3_x + c*pow(r,-1.0)*vis3_xy - p_xyz*pow(r,-1.0) + p_xz*pow(r,-2.0)*r_y + p_yz*pow(r,-2.0)*r_x - 2.0*p_z*pow(r,-3.0)*r_x*r_y + p_z*pow(r,-2.0)*r_xy - u*w_xxy - u_x*w_xy - u_xy*w_x - u_y*w_xx - v*w_xyy - v_x*w_yy - v_xy*w_y - v_y*w_xy - w*w_xyz - w_x*w_yz - w_xy*w_z - w_xz*w_y
p_txy = -c2*u*vis1_xy - c2*u_x*vis1_y - c2*u_xy*vis1 - c2*u_y*vis1_x - c2*v*vis2_xy - c2*v_x*vis2_y - c2*v_xy*vis2 - c2*v_y*vis2_x - c2*vis3*w_xy - c2*vis3_x*w_y - c2*vis3_xy*w - c2*vis3_y*w_x - gm*p*u_xxy - gm*p*v_xyy - gm*p*w_xyz - gm*p_x*u_xy - gm*p_x*v_yy - gm*p_x*w_yz - gm*p_xy*u_x - gm*p_xy*v_y - gm*p_xy*w_z - gm*p_y*u_xx - gm*p_y*v_xy - gm*p_y*w_xz - p_x*u_xy - p_xy*u_x - p_xy*v_y - p_xyz*w - p_xyy*v - p_xz*w_y - p_xx*u_y - p_xxy*u - p_y*v_xy - p_yz*w_x - p_yy*v_x - p_z*w_xy

r_tyz = -r*u_xyz - r*v_yyz - r*w_yzz - r_x*u_yz - r_xy*u_z - r_xyz*u - r_xz*u_y - r_y*u_xz - 2.0*r_y*v_yz - r_y*w_zz - r_yz*u_x - 2.0*r_yz*v_y - 2.0*r_yz*w_z - r_yzz*w - r_yy*v_z - r_yyz*v - r_z*u_xy - r_z*v_yy - 2.0*r_z*w_yz - r_zz*w_y
u_tyz = 2.0*c*pow(r,-3.0)*r_y*r_z*vis1 - c*pow(r,-2.0)*r_y*vis1_z - c*pow(r,-2.0)*r_yz*vis1 - c*pow(r,-2.0)*r_z*vis1_y + c*pow(r,-1.0)*vis1_yz - 2.0*p_x*pow(r,-3.0)*r_y*r_z + p_x*pow(r,-2.0)*r_yz + p_xy*pow(r,-2.0)*r_z - p_xyz*pow(r,-1.0) + p_xz*pow(r,-2.0)*r_y - u*u_xyz - u_x*u_yz - u_xy*u_z - u_xz*u_y - u_y*v_yz - u_yz*v_y - u_yz*w_z - u_yzz*w - u_yy*v_z - u_yyz*v - u_z*w_yz - u_zz*w_y
v_tyz = 2.0*c*pow(r,-3.0)*r_y*r_z*vis2 - c*pow(r,-2.0)*r_y*vis2_z - c*pow(r,-2.0)*r_yz*vis2 - c*pow(r,-2.0)*r_z*vis2_y + c*pow(r,-1.0)*vis2_yz - 2.0*p_y*pow(r,-3.0)*r_y*r_z + p_y*pow(r,-2.0)*r_yz + p_yz*pow(r,-2.0)*r_y + p_yy*pow(r,-2.0)*r_z - p_yyz*pow(r,-1.0) - u*v_xyz - u_y*v_xz - u_yz*v_x - u_z*v_xy - v*v_yyz - 2.0*v_y*v_yz - v_yz*w_z - v_yzz*w - v_yy*v_z - v_z*w_yz - v_zz*w_y
w_tyz = 2.0*c*pow(r,-3.0)*r_y*r_z*vis3 - c*pow(r,-2.0)*r_y*vis3_z - c*pow(r,-2.0)*r_yz*vis3 - c*pow(r,-2.0)*r_z*vis3_y + c*pow(r,-1.0)*vis3_yz + p_yz*pow(r,-2.0)*r_z - p_yzz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*r_y*r_z + p_z*pow(r,-2.0)*r_yz + p_zz*pow(r,-2.0)*r_y - u*w_xyz - u_y*w_xz - u_yz*w_x - u_z*w_xy - v*w_yyz - v_y*w_yz - v_yz*w_y - v_z*w_yy - w*w_yzz - w_y*w_zz - 2.0*w_yz*w_z
p_tyz = -c2*u*vis1_yz - c2*u_y*vis1_z - c2*u_yz*vis1 - c2*u_z*vis1_y - c2*v*vis2_yz - c2*v_y*vis2_z - c2*v_yz*vis2 - c2*v_z*vis2_y - c2*vis3*w_yz - c2*vis3_y*w_z - c2*vis3_yz*w - c2*vis3_z*w_y - gm*p*u_xyz - gm*p*v_yyz - gm*p*w_yzz - gm*p_y*u_xz - gm*p_y*v_yz - gm*p_y*w_zz - gm*p_yz*u_x - gm*p_yz*v_y - gm*p_yz*w_z - gm*p_z*u_xy - gm*p_z*v_yy - gm*p_z*w_yz - p_x*u_yz - p_xy*u_z - p_xyz*u - p_xz*u_y - p_y*v_yz - p_yz*v_y - p_yz*w_z - p_yzz*w - p_yy*v_z - p_yyz*v - p_z*w_yz - p_zz*w_y

r_txz = -r*u_xxz - r*v_xyz - r*w_xzz - 2.0*r_x*u_xz - r_x*v_yz - r_x*w_zz - r_xy*v_z - r_xyz*v - 2.0*r_xz*u_x - r_xz*v_y - 2.0*r_xz*w_z - r_xzz*w - r_xx*u_z - r_xxz*u - r_y*v_xz - r_yz*v_x - r_z*u_xx - r_z*v_xy - 2.0*r_z*w_xz - r_zz*w_x
u_txz = 2.0*c*pow(r,-3.0)*r_x*r_z*vis1 - c*pow(r,-2.0)*r_x*vis1_z - c*pow(r,-2.0)*r_xz*vis1 - c*pow(r,-2.0)*r_z*vis1_x + c*pow(r,-1.0)*vis1_xz - 2.0*p_x*pow(r,-3.0)*r_x*r_z + p_x*pow(r,-2.0)*r_xz + p_xz*pow(r,-2.0)*r_x + p_xx*pow(r,-2.0)*r_z - p_xxz*pow(r,-1.0) - u*u_xxz - 2.0*u_x*u_xz - u_xy*v_z - u_xyz*v - u_xz*w_z - u_xzz*w - u_xx*u_z - u_y*v_xz - u_yz*v_x - u_z*w_xz - u_zz*w_x
v_txz = 2.0*c*pow(r,-3.0)*r_x*r_z*vis2 - c*pow(r,-2.0)*r_x*vis2_z - c*pow(r,-2.0)*r_xz*vis2 - c*pow(r,-2.0)*r_z*vis2_x + c*pow(r,-1.0)*vis2_xz + p_xy*pow(r,-2.0)*r_z - p_xyz*pow(r,-1.0) - 2.0*p_y*pow(r,-3.0)*r_x*r_z + p_y*pow(r,-2.0)*r_xz + p_yz*pow(r,-2.0)*r_x - u*v_xxz - u_x*v_xz - u_xz*v_x - u_z*v_xx - v*v_xyz - v_x*v_yz - v_xy*v_z - v_xz*v_y - v_xz*w_z - v_xzz*w - v_z*w_xz - v_zz*w_x
w_txz = 2.0*c*pow(r,-3.0)*r_x*r_z*vis3 - c*pow(r,-2.0)*r_x*vis3_z - c*pow(r,-2.0)*r_xz*vis3 - c*pow(r,-2.0)*r_z*vis3_x + c*pow(r,-1.0)*vis3_xz + p_xz*pow(r,-2.0)*r_z - p_xzz*pow(r,-1.0) - 2.0*p_z*pow(r,-3.0)*r_x*r_z + p_z*pow(r,-2.0)*r_xz + p_zz*pow(r,-2.0)*r_x - u*w_xxz - u_x*w_xz - u_xz*w_x - u_z*w_xx - v*w_xyz - v_x*w_yz - v_xz*w_y - v_z*w_xy - w*w_xzz - w_x*w_zz - 2.0*w_xz*w_z
p_txz = -c2*u*vis1_xz - c2*u_x*vis1_z - c2*u_xz*vis1 - c2*u_z*vis1_x - c2*v*vis2_xz - c2*v_x*vis2_z - c2*v_xz*vis2 - c2*v_z*vis2_x - c2*vis3*w_xz - c2*vis3_x*w_z - c2*vis3_xz*w - c2*vis3_z*w_x - gm*p*u_xxz - gm*p*v_xyz - gm*p*w_xzz - gm*p_x*u_xz - gm*p_x*v_yz - gm*p_x*w_zz - gm*p_xz*u_x - gm*p_xz*v_y - gm*p_xz*w_z - gm*p_z*u_xx - gm*p_z*v_xy - gm*p_z*w_xz - p_x*u_xz - p_xy*v_z - p_xyz*v - p_xz*u_x - p_xz*w_z - p_xzz*w - p_xx*u_z - p_xxz*u - p_y*v_xz - p_yz*v_x - p_z*w_xz - p_zz*w_x

r_tt = -r*u_tx - r*v_ty - r*w_tz - r_t*u_x - r_t*v_y - r_t*w_z - r_tx*u - r_ty*v - r_tz*w - r_x*u_t - r_y*v_t - r_z*w_t
u_tt = -c*pow(r,-2.0)*r_t*vis1 + c*pow(r,-1.0)*vis1_t - p_tx*pow(r,-1.0) + p_x*pow(r,-2.0)*r_t - u*u_tx - u_t*u_x - u_ty*v - u_tz*w - u_y*v_t - u_z*w_t
v_tt = -c*pow(r,-2.0)*r_t*vis2 + c*pow(r,-1.0)*vis2_t - p_ty*pow(r,-1.0) + p_y*pow(r,-2.0)*r_t - u*v_tx - u_t*v_x - v*v_ty - v_t*v_y - v_tz*w - v_z*w_t
w_tt = -c*pow(r,-2.0)*r_t*vis3 + c*pow(r,-1.0)*vis3_t - p_tz*pow(r,-1.0) + p_z*pow(r,-2.0)*r_t - u*w_tx - u_t*w_x - v*w_ty - v_t*w_y - w*w_tz - w_t*w_z
p_tt = -c2*u*vis1_t - c2*u_t*vis1 - c2*v*vis2_t - c2*v_t*vis2 - c2*vis3*w_t - c2*vis3_t*w - gm*p*u_tx - gm*p*v_ty - gm*p*w_tz - gm*p_t*u_x - gm*p_t*v_y - gm*p_t*w_z - p_tx*u - p_ty*v - p_tz*w - p_x*u_t - p_y*v_t - p_z*w_t
