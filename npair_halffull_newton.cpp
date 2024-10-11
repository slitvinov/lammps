#include "npair_halffull_newton.h"
#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalffullNewton::NPairHalffullNewton(LAMMPS *lmp) : NPair(lmp) {}
void NPairHalffullNewton::build(NeighList *list)
{
  int i, j, ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;
  int inum = 0;
  ipage->reset();
  for (ii = 0; ii < inum_full; ii++) {
    n = 0;
    neighptr = ipage->vget();
    i = ilist_full[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];
    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j < nlocal) {
        if (i > j) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }
      neighptr[n++] = joriginal;
    }
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  list->inum = inum;
}
