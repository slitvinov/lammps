#include "npair_halffull_newtoff.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalffullNewtoff::NPairHalffullNewtoff(LAMMPS *lmp) : NPair(lmp) {}
void NPairHalffullNewtoff::build(NeighList *list)
{
  int i, j, ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
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
  for (ii = 0; ii < inum_full; ii++) {
    n = 0;
    neighptr = ipage->vget();
    i = ilist_full[ii];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];
    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j > i) neighptr[n++] = joriginal;
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
