#include "npair_half_size_bin_newton_tri.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalfSizeBinNewtonTri::NPairHalfSizeBinNewtonTri(LAMMPS *lmp) :
  NPair(lmp) {}
void NPairHalfSizeBinNewtonTri::build(NeighList *list)
{
  int i,j,jh,k,n,ibin,which;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;
  int history = list->history;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  int mask_history = 1 << HISTBITS;
  int inum = 0;
  ipage->reset();
  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    ibin = atom2bin[i];
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  list->inum = inum;
}
