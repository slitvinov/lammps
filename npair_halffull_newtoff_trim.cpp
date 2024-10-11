#include "npair_halffull_newtoff_trim.h"
#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalffullNewtoffTrim::NPairHalffullNewtoffTrim(LAMMPS *lmp) : NPair(lmp) {}
void NPairHalffullNewtoffTrim::build(NeighList *list)
{
  int i, j, ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq;
  double **x = atom->x;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;
  if (list->ghost) inum_full += list->listfull->gnum;
  int inum = 0;
  ipage->reset();
  double cutsq_custom = cutoff_custom * cutoff_custom;
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
      if (j > i) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq > cutsq_custom) continue;
        neighptr[n++] = joriginal;
      }
    }
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  list->inum = inum;
  if (list->ghost) list->gnum = list->listfull->gnum;
}
