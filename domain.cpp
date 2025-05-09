#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <cstdio>
#include <string>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "domain.h"
#include "lmptype.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "region_block.h"
#include "universe.h"
#include "update.h"
using namespace LAMMPS_NS;
#define BIG 1.0e20
#define SMALL 1.0e-4
#define DELTAREGION 4
#define BONDSTRETCH 1.1
Domain::Domain(LAMMPS *lmp) : Pointers(lmp) {
  box_exist = 0;
  box_change = 0;
  dimension = 3;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;
  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;
  minxlo = minxhi = 0.0;
  minylo = minyhi = 0.0;
  minzlo = minzhi = 0.0;
  triclinic = 0;
  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  xy = xz = yz = 0.0;
  h[3] = h[4] = h[5] = 0.0;
  h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
  h_rate[0] = h_rate[1] = h_rate[2] = h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;
  copymode = 0;
}
void Domain::init() {
  box_change_size = box_change_shape = box_change_domain = 0;
  int box_change_x = 0, box_change_y = 0, box_change_z = 0;
  int box_change_yz = 0, box_change_xz = 0, box_change_xy = 0;
  box_change = 0;
}
void Domain::set_initial_box(int expandflag) {
  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);
}
void Domain::set_global_box() {
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];
  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0 / h[0];
  h_inv[1] = 1.0 / h[1];
  h_inv[2] = 1.0 / h[2];
  prd_half[0] = xprd_half = 0.5 * xprd;
  prd_half[1] = yprd_half = 0.5 * yprd;
  prd_half[2] = zprd_half = 0.5 * zprd;
}
void Domain::set_local_box() {
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;
  sublo[0] = boxlo[0] + xprd * xsplit[myloc[0]];
  if (myloc[0] < procgrid[0] - 1)
    subhi[0] = boxlo[0] + xprd * xsplit[myloc[0] + 1];
  else
    subhi[0] = boxhi[0];
  sublo[1] = boxlo[1] + yprd * ysplit[myloc[1]];
  if (myloc[1] < procgrid[1] - 1)
    subhi[1] = boxlo[1] + yprd * ysplit[myloc[1] + 1];
  else
    subhi[1] = boxhi[1];
  sublo[2] = boxlo[2] + zprd * zsplit[myloc[2]];
  if (myloc[2] < procgrid[2] - 1)
    subhi[2] = boxlo[2] + zprd * zsplit[myloc[2] + 1];
  else
    subhi[2] = boxhi[2];
}
void Domain::reset_box() {
  if (atom->natoms == 0)
    return;
  set_global_box();
  set_local_box();
}
void Domain::pbc() {
  int nlocal = atom->nlocal;
  if (!nlocal)
    return;
  int i;
  imageint idim, otherdims;
  double *lo, *hi, *period;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double *coord;
  int n3 = 3 * nlocal;
  coord = &x[0][0];
  int flag = 0;
  for (i = 0; i < n3; i++)
    if (!std::isfinite(*coord++))
      flag = 1;
  lo = boxlo;
  hi = boxhi;
  period = prd;
}
