#include "math_extra.h"
#include <cstdio>
#include <cstring>
namespace MathExtra {
void write3(const double mat[3][3])
{
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) printf("%g ",mat[i][j]);
    printf("\n");
  }
}
int mldivide3(const double m[3][3], const double *v, double *ans)
{
  double aug[3][4];
  for (unsigned i = 0; i < 3; i++) {
    aug[i][3] = v[i];
    for (unsigned j = 0; j < 3; j++) aug[i][j] = m[i][j];
  }
  for (unsigned i = 0; i < 2; i++) {
    unsigned p = i;
    for (unsigned j = i+1; j < 3; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memmove(aug[i],aug[j],4*sizeof(double));
        memcpy(aug[j],tempv,4*sizeof(double));
      }
    }
    while (p < 3 && aug[p][i] == 0.0) p++;
    if (p == 3) return 1;
    else
      if (p != i) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memmove(aug[i],aug[p],4*sizeof(double));
        memcpy(aug[p],tempv,4*sizeof(double));
      }
    for (unsigned j = i+1; j < 3; j++) {
      double n = aug[j][i]/aug[i][i];
      for (unsigned k=i+1; k<4; k++) aug[j][k]-=n*aug[i][k];
    }
  }
  if (aug[2][2] == 0.0) return 1;
  ans[2] = aug[2][3]/aug[2][2];
  for (int i = 1; i >= 0; i--) {
    double sumax = 0.0;
    for (unsigned j = i+1; j < 3; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][3]-sumax) / aug[i][i];
  }
  return 0;
}
void richardson(double *q, double *m, double *w, double *moments, double dtq)
{
  double wq[4];
  MathExtra::vecquat(w,q,wq);
  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtra::qnormalize(qfull);
  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);
  MathExtra::mq_to_omega(m,qhalf,moments,w);
  MathExtra::vecquat(w,qhalf,wq);
  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);
  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];
  MathExtra::qnormalize(q);
}
void richardson_sphere(double *q, double *w, double dtq)
{
  double wq[4];
  MathExtra::vecquat(w,q,wq);
  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtra::qnormalize(qfull);
  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);
  MathExtra::vecquat(w,qhalf,wq);
  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);
  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];
  MathExtra::qnormalize(q);
}
void no_squish_rotate(int k, double *p, double *q, double *inertia,
                      double dt)
{
  double phi,c_phi,s_phi,kp[4],kq[4];
  if (k == 1) {
    kq[0] = -q[1]; kp[0] = -p[1];
    kq[1] = q[0]; kp[1] = p[0];
    kq[2] = q[3]; kp[2] = p[3];
    kq[3] = -q[2]; kp[3] = -p[2];
  } else if (k == 2) {
    kq[0] = -q[2]; kp[0] = -p[2];
    kq[1] = -q[3]; kp[1] = -p[3];
    kq[2] = q[0]; kp[2] = p[0];
    kq[3] = q[1]; kp[3] = p[1];
  } else if (k == 3) {
    kq[0] = -q[3]; kp[0] = -p[3];
    kq[1] = q[2]; kp[1] = p[2];
    kq[2] = -q[1]; kp[2] = -p[1];
    kq[3] = q[0]; kp[3] = p[0];
  }
  phi = p[0]*kq[0] + p[1]*kq[1] + p[2]*kq[2] + p[3]*kq[3];
  if (inertia[k-1] == 0.0) phi = 0.0;
  else phi /= 4.0 * inertia[k-1];
  c_phi = cos(dt * phi);
  s_phi = sin(dt * phi);
  p[0] = c_phi*p[0] + s_phi*kp[0];
  p[1] = c_phi*p[1] + s_phi*kp[1];
  p[2] = c_phi*p[2] + s_phi*kp[2];
  p[3] = c_phi*p[3] + s_phi*kp[3];
  q[0] = c_phi*q[0] + s_phi*kq[0];
  q[1] = c_phi*q[1] + s_phi*kq[1];
  q[2] = c_phi*q[2] + s_phi*kq[2];
  q[3] = c_phi*q[3] + s_phi*kq[3];
}
void angmom_to_omega(double *m, double *ex, double *ey, double *ez,
                     double *idiag, double *w)
{
  double wbody[3];
  if (idiag[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / idiag[0];
  if (idiag[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / idiag[1];
  if (idiag[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / idiag[2];
  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
}
void mq_to_omega(double *m, double *q, double *moments, double *w)
{
  double wbody[3];
  double rot[3][3];
  MathExtra::quat_to_mat(q,rot);
  MathExtra::transpose_matvec(rot,m,wbody);
  if (moments[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] /= moments[0];
  if (moments[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] /= moments[1];
  if (moments[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] /= moments[2];
  MathExtra::matvec(rot,wbody,w);
}
void omega_to_angmom(double *w, double *ex, double *ey, double *ez,
                     double *idiag, double *m)
{
  double mbody[3];
  mbody[0] = (w[0]*ex[0] + w[1]*ex[1] + w[2]*ex[2]) * idiag[0];
  mbody[1] = (w[0]*ey[0] + w[1]*ey[1] + w[2]*ey[2]) * idiag[1];
  mbody[2] = (w[0]*ez[0] + w[1]*ez[1] + w[2]*ez[2]) * idiag[2];
  m[0] = mbody[0]*ex[0] + mbody[1]*ey[0] + mbody[2]*ez[0];
  m[1] = mbody[0]*ex[1] + mbody[1]*ey[1] + mbody[2]*ez[1];
  m[2] = mbody[0]*ex[2] + mbody[1]*ey[2] + mbody[2]*ez[2];
}
void exyz_to_q(double *ex, double *ey, double *ez, double *q)
{
  double q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0);
  double q1sq = q0sq - 0.5 * (ey[1] + ez[2]);
  double q2sq = q0sq - 0.5 * (ex[0] + ez[2]);
  double q3sq = q0sq - 0.5 * (ex[0] + ey[1]);
  if (q0sq >= 0.25) {
    q[0] = sqrt(q0sq);
    q[1] = (ey[2] - ez[1]) / (4.0*q[0]);
    q[2] = (ez[0] - ex[2]) / (4.0*q[0]);
    q[3] = (ex[1] - ey[0]) / (4.0*q[0]);
  } else if (q1sq >= 0.25) {
    q[1] = sqrt(q1sq);
    q[0] = (ey[2] - ez[1]) / (4.0*q[1]);
    q[2] = (ey[0] + ex[1]) / (4.0*q[1]);
    q[3] = (ex[2] + ez[0]) / (4.0*q[1]);
  } else if (q2sq >= 0.25) {
    q[2] = sqrt(q2sq);
    q[0] = (ez[0] - ex[2]) / (4.0*q[2]);
    q[1] = (ey[0] + ex[1]) / (4.0*q[2]);
    q[3] = (ez[1] + ey[2]) / (4.0*q[2]);
  } else if (q3sq >= 0.25) {
    q[3] = sqrt(q3sq);
    q[0] = (ex[1] - ey[0]) / (4.0*q[3]);
    q[1] = (ez[0] + ex[2]) / (4.0*q[3]);
    q[2] = (ez[1] + ey[2]) / (4.0*q[3]);
  }
  qnormalize(q);
}
void q_to_exyz(double *q, double *ex, double *ey, double *ez)
{
  ex[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  ex[1] = 2.0 * (q[1]*q[2] + q[0]*q[3]);
  ex[2] = 2.0 * (q[1]*q[3] - q[0]*q[2]);
  ey[0] = 2.0 * (q[1]*q[2] - q[0]*q[3]);
  ey[1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  ey[2] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
  ez[0] = 2.0 * (q[1]*q[3] + q[0]*q[2]);
  ez[1] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
  ez[2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}
void quat_to_mat(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];
  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;
  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;
  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}
void quat_to_mat_trans(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];
  mat[0][0] = w2+i2-j2-k2;
  mat[1][0] = twoij-twokw;
  mat[2][0] = twojw+twoik;
  mat[0][1] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[2][1] = twojk-twoiw;
  mat[0][2] = twoik-twojw;
  mat[1][2] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}
void inertia_ellipsoid(double *radii, double *quat, double mass,
                       double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];
  double idiag[3];
  quat_to_mat(quat,p);
  quat_to_mat_trans(quat,ptrans);
  idiag[0] = 0.2*mass * (radii[1]*radii[1] + radii[2]*radii[2]);
  idiag[1] = 0.2*mass * (radii[0]*radii[0] + radii[2]*radii[2]);
  idiag[2] = 0.2*mass * (radii[0]*radii[0] + radii[1]*radii[1]);
  diag_times3(idiag,ptrans,itemp);
  times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}
void inertia_line(double length, double theta, double mass, double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];
  double q[4],idiag[3];
  q[0] = cos(0.5*theta);
  q[1] = q[2] = 0.0;
  q[3] = sin(0.5*theta);
  MathExtra::quat_to_mat(q,p);
  MathExtra::quat_to_mat_trans(q,ptrans);
  idiag[0] = 0.0;
  idiag[1] = 1.0/12.0 * mass * length*length;
  idiag[2] = 1.0/12.0 * mass * length*length;
  MathExtra::diag_times3(idiag,ptrans,itemp);
  MathExtra::times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}
void inertia_triangle(double *v0, double *v1, double *v2,
                      double mass, double *inertia)
{
  double s[3][3] = {{2.0, 1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 2.0}};
  double v[3][3],sv[3][3],vtsv[3][3];
  double vvv[3],v1mv0[3],v2mv0[3],normal[3];
  v[0][0] = v0[0]; v[0][1] = v0[1]; v[0][2] = v0[2];
  v[1][0] = v1[0]; v[1][1] = v1[1]; v[1][2] = v1[2];
  v[2][0] = v2[0]; v[2][1] = v2[1]; v[2][2] = v2[2];
  times3(s,v,sv);
  transpose_times3(v,sv,vtsv);
  double sum = lensq3(v0) + lensq3(v1) + lensq3(v2);
  vvv[0] = v0[0] + v1[0] + v2[0];
  vvv[1] = v0[1] + v1[1] + v2[1];
  vvv[2] = v0[2] + v1[2] + v2[2];
  sum += lensq3(vvv);
  sub3(v1,v0,v1mv0);
  sub3(v2,v0,v2mv0);
  cross3(v1mv0,v2mv0,normal);
  double a = len3(normal);
  double inv24 = mass/24.0;
  inertia[0] = inv24*a*(sum-vtsv[0][0]);
  inertia[1] = inv24*a*(sum-vtsv[1][1]);
  inertia[2] = inv24*a*(sum-vtsv[2][2]);
  inertia[3] = -inv24*a*vtsv[1][2];
  inertia[4] = -inv24*a*vtsv[0][2];
  inertia[5] = -inv24*a*vtsv[0][1];
}
void inertia_triangle(double *idiag, double *quat, double ,
                      double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];
  quat_to_mat(quat,p);
  quat_to_mat_trans(quat,ptrans);
  diag_times3(idiag,ptrans,itemp);
  times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}
void BuildRxMatrix(double R[3][3], const double angle)
{
  const double angleSq = angle * angle;
  const double cosAngle = (1.0 - angleSq * 0.25) / (1.0 + angleSq * 0.25);
  const double sinAngle = angle / (1.0 + angleSq * 0.25);
  R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
  R[1][0] = 0.0; R[1][1] = cosAngle; R[1][2] = -sinAngle;
  R[2][0] = 0.0; R[2][1] = sinAngle; R[2][2] = cosAngle;
}
void BuildRyMatrix(double R[3][3], const double angle)
{
  const double angleSq = angle * angle;
  const double cosAngle = (1.0 - angleSq * 0.25) / (1.0 + angleSq * 0.25);
  const double sinAngle = angle / (1.0 + angleSq * 0.25);
  R[0][0] = cosAngle; R[0][1] = 0.0; R[0][2] = sinAngle;
  R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
  R[2][0] = -sinAngle; R[2][1] = 0.0; R[2][2] = cosAngle;
}
void BuildRzMatrix(double R[3][3], const double angle)
{
  const double angleSq = angle * angle;
  const double cosAngle = (1.0 - angleSq * 0.25) / (1.0 + angleSq * 0.25);
  const double sinAngle = angle / (1.0 + angleSq * 0.25);
  R[0][0] = cosAngle; R[0][1] = -sinAngle; R[0][2] = 0.0;
  R[1][0] = sinAngle; R[1][1] = cosAngle; R[1][2] = 0.0;
  R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
}
void tribbox(double *h, double radius, double *dist)
{
  double lx = h[0];
  double ly = h[1];
  double lz = h[2];
  double yz = h[3];
  double xz = h[4];
  double xy = h[5];
  dist[0] = radius * sqrt(ly*ly*lz*lz + ly*ly*xz*xz - 2.0*ly*xy*xz*yz +
                          lz*lz*xy*xy + xy*xy*yz*yz) / (lx*ly*lz);
  dist[1] = radius * sqrt(lz*lz + yz*yz) / (ly*lz);
  dist[2] = radius / lz;
}
}
