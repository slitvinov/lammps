#include <map>
#include <set>
#include <vector>
#include "npair.h"
#include "atom.h"
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
}
void NPair::build_setup() {
  if (nb)
    copy_bin_info();
  last_build = update->ntimestep;
}
