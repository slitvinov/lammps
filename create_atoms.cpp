#include <set>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include <exception>
#include <unordered_set>
#include "pointers.h"
#include "create_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "random_park.h"
#include "region.h"
using namespace LAMMPS_NS;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::THIRD;
static constexpr double BIG = 1.0e30;
static constexpr double EPSILON = 1.0e-6;
static constexpr double LB_FACTOR = 1.1;
static constexpr double INV_P_CONST = 0.7548777;
static constexpr double INV_SQ_P_CONST = 0.5698403;
static constexpr int DEFAULT_MAXTRY = 1000;
enum { BOX, REGION, SINGLE, RANDOM, MESH };
enum { ATOM, MOLECULE };
enum { COUNT, INSERT, INSERT_SELECTED };
enum { NONE, RATIO, SUBSET };
enum { BISECTION, QUASIRANDOM };
static constexpr const char *mesh_name[] = {"recursive bisection",
                                            "quasi-random"};
CreateAtoms::CreateAtoms(LAMMPS *lmp) : Command(lmp), basistype(nullptr) {}
void CreateAtoms::command(int narg, char **arg) {
  int latsty = domain->lattice->style;
  ntype = utils::inumeric(FLERR, arg[0], false, lmp);
  const char *meshfile;
  int iarg;
  style = RANDOM;
  nrandom = utils::bnumeric(FLERR, arg[2], false, lmp);
  seed = utils::inumeric(FLERR, arg[3], false, lmp);
  region = domain->get_region_by_id(arg[4]);
  region->init();
  region->prematch();
  iarg = 5;
  int scaleflag = 1;
  remapflag = 0;
  mode = ATOM;
  int molseed;
  ranmol = nullptr;
  vstr = xstr = ystr = zstr = nullptr;
  quat_user = 0;
  quatone[0] = quatone[1] = quatone[2] = quatone[3] = 0.0;
  subsetflag = NONE;
  int subsetseed;
  maxtry = DEFAULT_MAXTRY;
  radscale = 1.0;
  mesh_style = BISECTION;
  mesh_density = 1.0;
  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++)
    basistype[i] = ntype;
  ranlatt = nullptr;
  if (subsetflag != NONE)
    ranlatt = new RanMars(lmp, subsetseed + comm->me);
  xone[0] *= domain->lattice->xlattice;
  xone[1] *= domain->lattice->ylattice;
  xone[2] *= domain->lattice->zlattice;
  overlap *= domain->lattice->xlattice;
  triclinic = domain->triclinic;
  double epsilon[3];
  epsilon[0] = domain->prd[0] * EPSILON;
  epsilon[1] = domain->prd[1] * EPSILON;
  epsilon[2] = domain->prd[2] * EPSILON;
  sublo[0] = domain->sublo[0];
  subhi[0] = domain->subhi[0];
  sublo[1] = domain->sublo[1];
  subhi[1] = domain->subhi[1];
  sublo[2] = domain->sublo[2];
  subhi[2] = domain->subhi[2];
  MPI_Barrier(world);
  atom->nghost = 0;
  atom->avec->clear_bonus();
  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;
  if (style == RANDOM)
    add_random();
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();
  delete ranmol;
  delete ranlatt;
  delete[] basistype;
  delete[] vstr;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  MPI_Barrier(world);
}
void CreateAtoms::add_random() {
  double xlo, ylo, zlo, xhi, yhi, zhi, zmid;
  double delx, dely, delz, distsq, odistsq;
  double lamda[3], *coord;
  double *boxlo, *boxhi;
  auto random = new RanPark(lmp, seed);
  for (int ii = 0; ii < 30; ii++)
    random->uniform();
  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  zmid = zlo + 0.5 * (zhi - zlo);
  if (region && region->bboxflag) {
    xlo = MAX(xlo, region->extent_xlo);
    xhi = MIN(xhi, region->extent_xhi);
    ylo = MAX(ylo, region->extent_ylo);
    yhi = MIN(yhi, region->extent_yhi);
    zlo = MAX(zlo, region->extent_zlo);
    zhi = MIN(zhi, region->extent_zhi);
  }
  int ntry, success;
  bigint ninsert = 0;
  for (bigint i = 0; i < nrandom; i++) {
    success = 0;
    ntry = 0;
    while (ntry < maxtry) {
      ntry++;
      xone[0] = xlo + random->uniform() * (xhi - xlo);
      xone[1] = ylo + random->uniform() * (yhi - ylo);
      xone[2] = zlo + random->uniform() * (zhi - zlo);
      if (region && (region->match(xone[0], xone[1], xone[2]) == 0))
        continue;
      coord = xone;
      success = 1;
      break;
    }
    if (!success)
      continue;
    ninsert++;
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
        coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      if (mode == ATOM)
        atom->avec->create_atom(ntype, xone);
    }
  }
  delete random;
}
