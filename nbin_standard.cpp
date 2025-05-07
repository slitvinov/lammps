#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include <cstdio>
#include <string>
#include <cmath>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "nbin.h"
#include "nbin_standard.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
using namespace LAMMPS_NS;
#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100
NBinStandard::NBinStandard(LAMMPS *lmp) : NBin(lmp) {}
void NBinStandard::bin_atoms_setup(int nall) {
  if (mbins > maxbin) {
    maxbin = mbins;
    memory->destroy(binhead);
    memory->create(binhead, maxbin, "neigh:binhead");
  }
  if (nall > maxatom) {
    maxatom = nall;
    memory->destroy(bins);
    memory->create(bins, maxatom, "neigh:bins");
    memory->destroy(atom2bin);
    memory->create(atom2bin, maxatom, "neigh:atom2bin");
  }
}
void NBinStandard::setup_bins(int style) {
  double bbox[3], bsubboxlo[3], bsubboxhi[3];
  double *cutghost = comm->cutghost;
  bsubboxlo[0] = domain->sublo[0] - cutghost[0];
  bsubboxlo[1] = domain->sublo[1] - cutghost[1];
  bsubboxlo[2] = domain->sublo[2] - cutghost[2];
  bsubboxhi[0] = domain->subhi[0] + cutghost[0];
  bsubboxhi[1] = domain->subhi[1] + cutghost[1];
  bsubboxhi[2] = domain->subhi[2] + cutghost[2];
  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];
  double binsize_optimal;
  binsize_optimal = binsize_user;
  if (binsize_optimal == 0.0)
    binsize_optimal = bbox[0];
  double binsizeinv = 1.0 / binsize_optimal;
  nbinx = static_cast<int>(bbox[0] * binsizeinv);
  nbiny = static_cast<int>(bbox[1] * binsizeinv);
  nbinz = static_cast<int>(bbox[2] * binsizeinv);
  binsizex = bbox[0] / nbinx;
  binsizey = bbox[1] / nbiny;
  binsizez = bbox[2] / nbinz;
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;
  int mbinxhi, mbinyhi, mbinzhi;
  double coord;
  coord = bsubboxlo[0] - SMALL * bbox[0];
  mbinxlo = static_cast<int>((coord - bboxlo[0]) * bininvx);
  if (coord < bboxlo[0])
    mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL * bbox[0];
  mbinxhi = static_cast<int>((coord - bboxlo[0]) * bininvx);
  coord = bsubboxlo[1] - SMALL * bbox[1];
  mbinylo = static_cast<int>((coord - bboxlo[1]) * bininvy);
  if (coord < bboxlo[1])
    mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL * bbox[1];
  mbinyhi = static_cast<int>((coord - bboxlo[1]) * bininvy);
  if (dimension == 3) {
    coord = bsubboxlo[2] - SMALL * bbox[2];
    mbinzlo = static_cast<int>((coord - bboxlo[2]) * bininvz);
    if (coord < bboxlo[2])
      mbinzlo = mbinzlo - 1;
    coord = bsubboxhi[2] + SMALL * bbox[2];
    mbinzhi = static_cast<int>((coord - bboxlo[2]) * bininvz);
  }
  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;
  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;
  if (dimension == 3) {
    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;
  } else
    mbinzlo = mbinzhi = 0;
  mbinz = mbinzhi - mbinzlo + 1;
  bigint bbin = ((bigint)mbinx) * ((bigint)mbiny) * ((bigint)mbinz) + 1;
  mbins = bbin;
}
void NBinStandard::bin_atoms() {
  int i, ibin;
  last_bin = update->ntimestep;
  for (i = 0; i < mbins; i++)
    binhead[i] = -1;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  for (i = nall - 1; i >= 0; i--) {
    ibin = coord2bin(x[i]);
    atom2bin[i] = ibin;
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }
}
