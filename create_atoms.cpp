#include <set>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include <exception>
#include <unordered_set>
#include <string>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "command.h"
#include "create_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
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
CreateAtoms::CreateAtoms(LAMMPS *lmp) : Command(lmp) {}
void CreateAtoms::command(int narg, char **arg) {
  ntype = utils::inumeric(FLERR, arg[0], false, lmp);
  int iarg;
  nrandom = utils::bnumeric(FLERR, arg[2], false, lmp);
  seed = utils::inumeric(FLERR, arg[3], false, lmp);
  region = domain->get_region_by_id(arg[4]);
  region->init();
  iarg = 5;
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
  add_random();
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();
  MPI_Barrier(world);
}
void CreateAtoms::add_random() {
  double xlo, ylo, zlo, xhi, yhi, zhi;
  double *coord;
  double *boxlo, *boxhi, xone[3];
  RanPark *random = new RanPark(lmp, seed);
  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  for (bigint i = 0; i < nrandom; i++) {
    xone[0] = xlo + random->uniform() * (xhi - xlo);
    xone[1] = ylo + random->uniform() * (yhi - ylo);
    xone[2] = zlo + random->uniform() * (zhi - zlo);
    coord = xone;
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
        coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      atom->avec->create_atom(ntype, xone);
    }
  }
  delete random;
}
