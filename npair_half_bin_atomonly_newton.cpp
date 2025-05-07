#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <mpi.h>
#include <string>
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "npair.h"
#include "npair_half_bin_atomonly_newton.h"
#include "atom.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalfBinAtomonlyNewton::NPairHalfBinAtomonlyNewton(LAMMPS *lmp)
    : NPair(lmp) {}
void NPairHalfBinAtomonlyNewton::build(NeighList *list) {
  int i, j, k, n, itype, jtype, ibin;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neighptr;
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (includegroup)
    nlocal = atom->nfirst;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  int inum = 0;
  ipage->reset();
  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
        if (x[j][2] < ztmp)
          continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp)
            continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp)
            continue;
        }
      }
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq <= cutneighsq[itype][jtype])
        neighptr[n++] = j;
    }
    ibin = atom2bin[i];
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
  }
  list->inum = inum;
}
