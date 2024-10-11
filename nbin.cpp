#include "nbin.h"
#include <cmath>
#include "neighbor.h"
#include "neigh_request.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
using namespace LAMMPS_NS;
NBin::NBin(LAMMPS *lmp) : Pointers(lmp)
{
  last_bin = -1;
  mbins = maxbin = maxatom = 0;
  binhead = nullptr;
  bins = nullptr;
  atom2bin = nullptr;
  nbinx_multi = nullptr; nbiny_multi = nullptr; nbinz_multi = nullptr;
  mbins_multi = nullptr;
  mbinx_multi = nullptr; mbiny_multi = nullptr, mbinz_multi = nullptr;
  mbinxlo_multi = nullptr;
  mbinylo_multi = nullptr;
  mbinzlo_multi = nullptr;
  binsizex_multi = nullptr; binsizey_multi = nullptr; binsizez_multi = nullptr;
  bininvx_multi = nullptr; bininvy_multi = nullptr; bininvz_multi = nullptr;
  binhead_multi = nullptr;
  maxbins_multi = nullptr;
  maxcollections = 0;
  neighbor->last_setup_bins = -1;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
}
NBin::~NBin()
{
  memory->destroy(binhead);
  memory->destroy(bins);
  memory->destroy(atom2bin);
  if (!binhead_multi) return;
  memory->destroy(nbinx_multi);
  memory->destroy(nbiny_multi);
  memory->destroy(nbinz_multi);
  memory->destroy(mbins_multi);
  memory->destroy(mbinx_multi);
  memory->destroy(mbiny_multi);
  memory->destroy(mbinz_multi);
  memory->destroy(mbinxlo_multi);
  memory->destroy(mbinylo_multi);
  memory->destroy(mbinzlo_multi);
  memory->destroy(binsizex_multi);
  memory->destroy(binsizey_multi);
  memory->destroy(binsizez_multi);
  memory->destroy(bininvx_multi);
  memory->destroy(bininvy_multi);
  memory->destroy(bininvz_multi);
  for (int n = 0; n < maxcollections; n++) {
    memory->destroy(binhead_multi[n]);
  }
  delete [] binhead_multi;
  memory->destroy(maxbins_multi);
}
void NBin::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}
void NBin::copy_neighbor_info()
{
  includegroup = neighbor->includegroup;
  cutneighmin = neighbor->cutneighmin;
  cutneighmax = neighbor->cutneighmax;
  binsizeflag = neighbor->binsizeflag;
  binsize_user = neighbor->binsize_user;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;
  ncollections = neighbor->ncollections;
  cutcollectionsq = neighbor->cutcollectionsq;
  if (cutoff_custom > 0.0) cutneighmax = cutoff_custom;
}
int NBin::coord2bin(double *x)
{
  int ix,iy,iz;
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");
  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;
  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;
  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;
  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}
int NBin::coord2bin_multi(double *x, int ic)
{
  int ix,iy,iz;
  int ibin;
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");
  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx_multi[ic]) + nbinx_multi[ic];
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_multi[ic]);
    ix = MIN(ix,nbinx_multi[ic]-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_multi[ic]) - 1;
  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy_multi[ic]) + nbiny_multi[ic];
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_multi[ic]);
    iy = MIN(iy,nbiny_multi[ic]-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_multi[ic]) - 1;
  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz_multi[ic]) + nbinz_multi[ic];
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_multi[ic]);
    iz = MIN(iz,nbinz_multi[ic]-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_multi[ic]) - 1;
  ibin = (iz-mbinzlo_multi[ic])*mbiny_multi[ic]*mbinx_multi[ic]
       + (iy-mbinylo_multi[ic])*mbinx_multi[ic]
       + (ix-mbinxlo_multi[ic]);
  return ibin;
}
