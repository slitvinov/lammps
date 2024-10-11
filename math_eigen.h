#ifndef LMP_MATH_EIGEN_H
#define LMP_MATH_EIGEN_H 
namespace MathEigen {
int jacobi3(double const *const *mat, double *eval, double **evec);
int jacobi3(double const mat[3][3], double *eval, double evec[3][3]);
}
#endif
