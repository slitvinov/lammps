#include "procmap.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"
#include "universe.h"
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
using namespace LAMMPS_NS;
#define MAXLINE 128
ProcMap::ProcMap(LAMMPS *lmp) : Pointers(lmp) {}
void ProcMap::onelevel_grid(int nprocs, int *user_procgrid, int *procgrid,
                            int otherflag, int other_style, int *other_procgrid,
                            int *other_coregrid) {
  int **factors;
  int npossible = factor(nprocs, nullptr);
  memory->create(factors, npossible, 3, "procmap:factors");
  npossible = factor(nprocs, factors);
  if (domain->dimension == 2)
    npossible = cull_2d(npossible, factors, 3);
  npossible = cull_user(npossible, factors, 3, user_procgrid);
  if (otherflag)
    npossible = cull_other(npossible, factors, 3, other_style, other_procgrid,
                           other_coregrid);
  if (npossible == 0)
    error->all(FLERR, "Could not create 3d grid of processors");
  best_factors(npossible, factors, procgrid, 1, 1, 1);
  memory->destroy(factors);
}
void ProcMap::cart_map(int reorder, int *procgrid, int *myloc,
                       int procneigh[3][2], int ***grid2proc) {
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;
  MPI_Cart_create(world, 3, procgrid, periods, reorder, &cartesian);
  MPI_Cart_get(cartesian, 3, procgrid, periods, myloc);
  MPI_Cart_shift(cartesian, 0, 1, &procneigh[0][0], &procneigh[0][1]);
  MPI_Cart_shift(cartesian, 1, 1, &procneigh[1][0], &procneigh[1][1]);
  MPI_Cart_shift(cartesian, 2, 1, &procneigh[2][0], &procneigh[2][1]);
  int coords[3];
  int i, j, k;
  for (i = 0; i < procgrid[0]; i++)
    for (j = 0; j < procgrid[1]; j++)
      for (k = 0; k < procgrid[2]; k++) {
        coords[0] = i;
        coords[1] = j;
        coords[2] = k;
        MPI_Cart_rank(cartesian, coords, &grid2proc[i][j][k]);
      }
  MPI_Comm_free(&cartesian);
}
int ProcMap::factor(int n, int **factors) {
  int i, j, nyz;
  int m = 0;
  for (i = 1; i <= n; i++) {
    if (n % i)
      continue;
    nyz = n / i;
    for (j = 1; j <= nyz; j++) {
      if (nyz % j)
        continue;
      if (factors) {
        factors[m][0] = i;
        factors[m][1] = j;
        factors[m][2] = nyz / j;
      }
      m++;
    }
  }
  return m;
}
int ProcMap::combine_factors(int n1, int **factors1, int n2, int **factors2,
                             int **factors) {
  int m = 0;
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      factors[m][0] = factors1[i][0] * factors2[j][0];
      factors[m][1] = factors1[i][1] * factors2[j][1];
      factors[m][2] = factors1[i][2] * factors2[j][2];
      factors[m][3] = j;
      m++;
    }
  return n1 * n2;
}
int ProcMap::cull_2d(int n, int **factors, int m) {
  int i = 0;
  while (i < n) {
    if (factors[i][2] != 1) {
      for (int j = 0; j < m; j++)
        factors[i][j] = factors[n - 1][j];
      n--;
    } else
      i++;
  }
  return n;
}
int ProcMap::cull_user(int n, int **factors, int m, int *user_factors) {
  int i = 0;
  while (i < n) {
    int flag = 0;
    if (user_factors[0] && factors[i][0] != user_factors[0])
      flag = 1;
    if (user_factors[1] && factors[i][1] != user_factors[1])
      flag = 1;
    if (user_factors[2] && factors[i][2] != user_factors[2])
      flag = 1;
    if (flag) {
      for (int j = 0; j < m; j++)
        factors[i][j] = factors[n - 1][j];
      n--;
    } else
      i++;
  }
  return n;
}
int ProcMap::cull_other(int n, int **factors, int m, int other_style,
                        int *other_procgrid, int *other_coregrid) {
  int i = 0;
  while (i < n) {
    if (other_style == Comm::MULTIPLE) {
      int flag = 0;
      if ((other_procgrid[0] / other_coregrid[0]) % factors[i][0])
        flag = 1;
      if ((other_procgrid[1] / other_coregrid[1]) % factors[i][1])
        flag = 1;
      if ((other_procgrid[2] / other_coregrid[2]) % factors[i][2])
        flag = 1;
      if (flag) {
        for (int j = 0; j < m; j++)
          factors[i][j] = factors[n - 1][j];
        n--;
      } else
        i++;
    }
  }
  return n;
}
int ProcMap::best_factors(int npossible, int **factors, int *best, const int sx,
                          const int sy, const int sz) {
  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx * sy);
    area[1] = domain->xprd * domain->zprd / (sx * sz);
    area[2] = domain->yprd * domain->zprd / (sy * sz);
  } else {
    double *h = domain->h;
    double a[3], b[3], c[3];
    a[0] = h[0];
    a[1] = 0.0;
    a[2] = 0.0;
    b[0] = h[5];
    b[1] = h[1];
    b[2] = 0.0;
    MathExtra::cross3(a, b, c);
    area[0] = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / (sx * sy);
    a[0] = h[0];
    a[1] = 0.0;
    a[2] = 0.0;
    b[0] = h[4];
    b[1] = h[3];
    b[2] = h[2];
    MathExtra::cross3(a, b, c);
    area[1] = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / (sx * sz);
    a[0] = h[5];
    a[1] = h[1];
    a[2] = 0.0;
    b[0] = h[4];
    b[1] = h[3];
    b[2] = h[2];
    MathExtra::cross3(a, b, c);
    area[2] = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / (sy * sz);
  }
  int index;
  double surf;
  double bestsurf = 2.0 * (area[0] + area[1] + area[2]);
  for (int m = 0; m < npossible; m++) {
    surf = area[0] / factors[m][0] / factors[m][1] +
           area[1] / factors[m][0] / factors[m][2] +
           area[2] / factors[m][1] / factors[m][2];
    if (surf < bestsurf) {
      bestsurf = surf;
      best[0] = factors[m][0];
      best[1] = factors[m][1];
      best[2] = factors[m][2];
      index = m;
    }
  }
  return index;
}
