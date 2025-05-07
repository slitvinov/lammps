#include <unordered_set>
#include <map>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
#include <cstdio>
#include <mpi.h>
#include "pointers.h"
#include "procmap.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "universe.h"
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
  npossible = cull_user(npossible, factors, 3, user_procgrid);
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
int ProcMap::best_factors(int npossible, int **factors, int *best, const int sx,
                          const int sy, const int sz) {
  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx * sy);
    area[1] = domain->xprd * domain->zprd / (sx * sz);
    area[2] = domain->yprd * domain->zprd / (sy * sz);
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
