#include "npair.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "nbin.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include <cmath>
using namespace LAMMPS_NS;
NPair::NPair(LAMMPS *lmp) : Pointers(lmp), nb(nullptr), bins(nullptr) {
  last_build = -1;
  mycutneighsq = nullptr;
  copymode = 0;
  execution_space = Host;
}
NPair::~NPair() {
  if (copymode)
    return;
  memory->destroy(mycutneighsq);
}
void NPair::post_constructor(NeighRequest *nrq) {
  cutoff_custom = 0.0;
  if (nrq->cut)
    cutoff_custom = nrq->cutoff;
}
void NPair::copy_neighbor_info() {
  includegroup = neighbor->includegroup;
  exclude = neighbor->exclude;
  skin = neighbor->skin;
  cutneighsq = neighbor->cutneighsq;
  cutneighghostsq = neighbor->cutneighghostsq;
  cut_inner_sq = neighbor->cut_inner_sq;
  cut_middle_sq = neighbor->cut_middle_sq;
  cut_middle_inside_sq = neighbor->cut_middle_inside_sq;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;
  nex_type = neighbor->nex_type;
  ex1_type = neighbor->ex1_type;
  ex2_type = neighbor->ex2_type;
  ex_type = neighbor->ex_type;
  nex_group = neighbor->nex_group;
  ex1_group = neighbor->ex1_group;
  ex2_group = neighbor->ex2_group;
  ex1_bit = neighbor->ex1_bit;
  ex2_bit = neighbor->ex2_bit;
  nex_mol = neighbor->nex_mol;
  ex_mol_group = neighbor->ex_mol_group;
  ex_mol_bit = neighbor->ex_mol_bit;
  ex_mol_intra = neighbor->ex_mol_intra;
  special_flag = neighbor->special_flag;
  ncollections = neighbor->ncollections;
  cutcollectionsq = neighbor->cutcollectionsq;
  if (cutoff_custom > 0.0) {
    memory->destroy(mycutneighsq);
    int n = atom->ntypes;
    memory->create(mycutneighsq, n + 1, n + 1, "npair:cutneighsq");
    int i, j;
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        mycutneighsq[i][j] = cutoff_custom * cutoff_custom;
    cutneighsq = mycutneighsq;
  }
}
void NPair::copy_bin_info() {
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;
  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;
  atom2bin = nb->atom2bin;
  bins = nb->bins;
  binhead = nb->binhead;
  nbinx_multi = nb->nbinx_multi;
  nbiny_multi = nb->nbiny_multi;
  nbinz_multi = nb->nbinz_multi;
  mbins_multi = nb->mbins_multi;
  mbinx_multi = nb->mbinx_multi;
  mbiny_multi = nb->mbiny_multi;
  mbinz_multi = nb->mbinz_multi;
  mbinxlo_multi = nb->mbinxlo_multi;
  mbinylo_multi = nb->mbinylo_multi;
  mbinzlo_multi = nb->mbinzlo_multi;
  bininvx_multi = nb->bininvx_multi;
  bininvy_multi = nb->bininvy_multi;
  bininvz_multi = nb->bininvz_multi;
  binhead_multi = nb->binhead_multi;
}
void NPair::build_setup() {
  if (nb)
    copy_bin_info();
  last_build = update->ntimestep;
}
int NPair::exclusion(int i, int j, int itype, int jtype, int *mask) const {
  int m;
  if (nex_type && ex_type[itype][jtype])
    return 1;
  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask[i] & ex1_bit[m] && mask[j] & ex2_bit[m])
        return 1;
      if (mask[i] & ex2_bit[m] && mask[j] & ex1_bit[m])
        return 1;
    }
  }
  return 0;
}
int NPair::coord2bin(double *x, int &ix, int &iy, int &iz) {
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR, "Non-numeric positions - simulation unstable");
  if (x[0] >= bboxhi[0])
    ix = static_cast<int>((x[0] - bboxhi[0]) * bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int>((x[0] - bboxlo[0]) * bininvx);
    ix = MIN(ix, nbinx - 1);
  } else
    ix = static_cast<int>((x[0] - bboxlo[0]) * bininvx) - 1;
  if (x[1] >= bboxhi[1])
    iy = static_cast<int>((x[1] - bboxhi[1]) * bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int>((x[1] - bboxlo[1]) * bininvy);
    iy = MIN(iy, nbiny - 1);
  } else
    iy = static_cast<int>((x[1] - bboxlo[1]) * bininvy) - 1;
  if (x[2] >= bboxhi[2])
    iz = static_cast<int>((x[2] - bboxhi[2]) * bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int>((x[2] - bboxlo[2]) * bininvz);
    iz = MIN(iz, nbinz - 1);
  } else
    iz = static_cast<int>((x[2] - bboxlo[2]) * bininvz) - 1;
  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  return iz * mbiny * mbinx + iy * mbinx + ix;
}
int NPair::coord2bin(double *x, int ic) {
  int ix, iy, iz;
  int ibin;
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR, "Non-numeric positions - simulation unstable");
  if (x[0] >= bboxhi[0])
    ix = static_cast<int>((x[0] - bboxhi[0]) * bininvx_multi[ic]) +
         nbinx_multi[ic];
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int>((x[0] - bboxlo[0]) * bininvx_multi[ic]);
    ix = MIN(ix, nbinx_multi[ic] - 1);
  } else
    ix = static_cast<int>((x[0] - bboxlo[0]) * bininvx_multi[ic]) - 1;
  if (x[1] >= bboxhi[1])
    iy = static_cast<int>((x[1] - bboxhi[1]) * bininvy_multi[ic]) +
         nbiny_multi[ic];
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int>((x[1] - bboxlo[1]) * bininvy_multi[ic]);
    iy = MIN(iy, nbiny_multi[ic] - 1);
  } else
    iy = static_cast<int>((x[1] - bboxlo[1]) * bininvy_multi[ic]) - 1;
  if (x[2] >= bboxhi[2])
    iz = static_cast<int>((x[2] - bboxhi[2]) * bininvz_multi[ic]) +
         nbinz_multi[ic];
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int>((x[2] - bboxlo[2]) * bininvz_multi[ic]);
    iz = MIN(iz, nbinz_multi[ic] - 1);
  } else
    iz = static_cast<int>((x[2] - bboxlo[2]) * bininvz_multi[ic]) - 1;
  ix -= mbinxlo_multi[ic];
  iy -= mbinylo_multi[ic];
  iz -= mbinzlo_multi[ic];
  ibin = iz * mbiny_multi[ic] * mbinx_multi[ic] + iy * mbinx_multi[ic] + ix;
  return ibin;
}
