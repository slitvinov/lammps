#include "npair_half_size_bin_newton.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairHalfSizeBinNewton::NPairHalfSizeBinNewton(LAMMPS *lmp) : NPair(lmp) {}
void NPairHalfSizeBinNewton::build(NeighList *list)
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
    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }
      if (exclude && exclusion(i,j,type[i],type[j],mask)) continue;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);
      if (rsq <= cutsq) {
        jh = j;
        if (history && rsq < radsum*radsum)
            jh = jh ^ mask_history;
 neighptr[n++] = jh;
      }
    }
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
