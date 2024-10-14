#include "math_eigen.h"
#include "math_eigen_impl.h"
#include <array>
#include <utility>
#include <vector>
using std::array;
using std::vector;
using namespace MathEigen;
typedef Jacobi<double, double *, double (*)[3], double const (*)[3]> Jacobi_v1;
typedef Jacobi<double, double *, double **, double const *const *> Jacobi_v2;
int MathEigen::jacobi3(double const mat[3][3], double *eval,
                       double evec[3][3]) {
  double mat_cpy[3][3] = {{mat[0][0], mat[0][1], mat[0][2]},
                          {mat[1][0], mat[1][1], mat[1][2]},
                          {mat[2][0], mat[2][1], mat[2][2]}};
  double *M[3] = {&(mat_cpy[0][0]), &(mat_cpy[1][0]), &(mat_cpy[2][0])};
  int midx[3];
  Jacobi_v1 ecalc3(3, M, midx);
  int ierror =
      ecalc3.Diagonalize(mat, eval, evec, Jacobi_v1::SORT_DECREASING_EVALS);
  for (int i = 0; i < 3; i++)
    for (int j = i + 1; j < 3; j++)
      std::swap(evec[i][j], evec[j][i]);
  return ierror;
}
int MathEigen::jacobi3(double const *const *mat, double *eval, double **evec) {
  double mat_cpy[3][3] = {{mat[0][0], mat[0][1], mat[0][2]},
                          {mat[1][0], mat[1][1], mat[1][2]},
                          {mat[2][0], mat[2][1], mat[2][2]}};
  double *M[3] = {&(mat_cpy[0][0]), &(mat_cpy[1][0]), &(mat_cpy[2][0])};
  int midx[3];
  Jacobi_v2 ecalc3(3, M, midx);
  int ierror =
      ecalc3.Diagonalize(mat, eval, evec, Jacobi_v2::SORT_DECREASING_EVALS);
  for (int i = 0; i < 3; i++)
    for (int j = i + 1; j < 3; j++)
      std::swap(evec[i][j], evec[j][i]);
  return ierror;
}
