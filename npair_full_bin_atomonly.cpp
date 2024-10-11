#include "npair_full_bin_atomonly.h"
#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairFullBinAtomonly::NPairFullBinAtomonly(LAMMPS *lmp) : NPair(lmp) {}
void NPairFullBinAtomonly::build(NeighList *list)
{
  int i, j, k, n, itype, jtype, ibin;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neighptr;
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = NULL;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;
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
    ibin = atom2bin[i];
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  list->inum = inum;
  list->gnum = 0;
}
