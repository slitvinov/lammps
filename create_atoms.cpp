#include "create_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "random_park.h"
#include "region.h"
#include <cmath>
#include <cstring>
#include <exception>
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
static constexpr const char *mesh_name[] = {"recursive bisection", "quasi-random"};
CreateAtoms::CreateAtoms(LAMMPS *lmp) : Command(lmp), basistype(nullptr) {}
void CreateAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Create_atoms command before simulation box is defined");
  if (modify->nfix_restart_peratom)
    error->all(FLERR, "Cannot create_atoms after reading restart file with per-atom info");
  int latsty = domain->lattice->style;
  if (latsty == Lattice::SQ || latsty == Lattice::SQ2 || latsty == Lattice::HEX)
    error->all(FLERR, "Lattice style incompatible with simulation dimension");
  if (narg < 2) utils::missing_cmd_args(FLERR, "create_atoms", error);
  ntype = utils::inumeric(FLERR, arg[0], false, lmp);
  const char *meshfile;
  int iarg;
  if (strcmp(arg[1], "random") == 0) {
    style = RANDOM;
    if (narg < 5) utils::missing_cmd_args(FLERR, "create_atoms random", error);
    nrandom = utils::bnumeric(FLERR, arg[2], false, lmp);
    if (nrandom < 0) error->all(FLERR, "Illegal create_atoms number of random atoms {}", nrandom);
    seed = utils::inumeric(FLERR, arg[3], false, lmp);
    if (seed <= 0) error->all(FLERR, "Illegal create_atoms random seed {}", seed);
    if (strcmp(arg[4], "NULL") == 0)
      region = nullptr;
    else {
      region = domain->get_region_by_id(arg[4]);
      if (!region) error->all(FLERR, "Create_atoms region {} does not exist", arg[4]);
      region->init();
      region->prematch();
    }
    iarg = 5;
  } else
    error->all(FLERR, "Unknown create_atoms command option {}", arg[1]);
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
  overlapflag = 0;
  maxtry = DEFAULT_MAXTRY;
  radscale = 1.0;
  radthresh = domain->lattice->xlattice;
  mesh_style = BISECTION;
  mesh_density = 1.0;
  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) basistype[i] = ntype;
  if (mode == ATOM) {
    if ((ntype <= 0) || (ntype > atom->ntypes))
      error->all(FLERR, "Invalid atom type in create_atoms command");
  }
  ranlatt = nullptr;
  if (subsetflag != NONE) ranlatt = new RanMars(lmp, subsetseed + comm->me);
  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR, "Incomplete use of variables in create_atoms command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR, "Incomplete use of variables in create_atoms command");
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
  double time1 = platform::walltime();
  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;
  if (style == RANDOM)
    add_random();
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT) error->all(FLERR, "Too many total atoms");
  if (atom->tag_enable) atom->tag_extend();
  atom->tag_check();
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
  delete ranmol;
  delete ranlatt;
  delete[] basistype;
  delete[] vstr;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  MPI_Barrier(world);
  if (comm->me == 0) {
    utils::logmesg(lmp, "Created {} atoms\n", atom->natoms - natoms_previous);
    if (scaleflag)
      domain->print_box("  using lattice units in ");
    else
      domain->print_box("  using box units in ");
    utils::logmesg(lmp, "  create_atoms CPU = {:.3f} seconds\n", platform::walltime() - time1);
  }
}
void CreateAtoms::add_random()
{
  double xlo, ylo, zlo, xhi, yhi, zhi, zmid;
  double delx, dely, delz, distsq, odistsq;
  double lamda[3], *coord;
  double *boxlo, *boxhi;
  if (overlapflag) {
    double odist = overlap;
    odistsq = odist * odist;
  }
  auto random = new RanPark(lmp, seed);
  for (int ii = 0; ii < 30; ii++) random->uniform();
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
  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR, "No overlap of box and region for create_atoms");
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
      if (domain->dimension == 2) xone[2] = zmid;
      if (region && (region->match(xone[0], xone[1], xone[2]) == 0)) continue;
      coord = xone;
      if (overlapflag) {
        double **x = atom->x;
        int nlocal = atom->nlocal;
        int reject = 0;
        for (int i = 0; i < nlocal; i++) {
          delx = xone[0] - x[i][0];
          dely = xone[1] - x[i][1];
          delz = xone[2] - x[i][2];
          domain->minimum_image(delx, dely, delz);
          distsq = delx * delx + dely * dely + delz * delz;
          if (distsq < odistsq) {
            reject = 1;
            break;
          }
        }
        int reject_any;
        MPI_Allreduce(&reject, &reject_any, 1, MPI_INT, MPI_MAX, world);
        if (reject_any) continue;
      }
      success = 1;
      break;
    }
    if (!success) continue;
    ninsert++;
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
        coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      if (mode == ATOM)
        atom->avec->create_atom(ntype, xone);
    }
  }
  if (ninsert < nrandom && comm->me == 0)
    error->warning(FLERR, "Only inserted {} particles out of {}", ninsert, nrandom);
  delete random;
}
