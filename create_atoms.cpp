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
#include "text_file_reader.h"
#include "variable.h"
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
  if (domain->dimension == 2) {
    if (latsty == Lattice::SC || latsty == Lattice::BCC || latsty == Lattice::FCC ||
        latsty == Lattice::HCP || latsty == Lattice::DIAMOND)
      error->all(FLERR, "Lattice style incompatible with simulation dimension");
  } else {
    if (latsty == Lattice::SQ || latsty == Lattice::SQ2 || latsty == Lattice::HEX)
      error->all(FLERR, "Lattice style incompatible with simulation dimension");
  }
  if (narg < 2) utils::missing_cmd_args(FLERR, "create_atoms", error);
  ntype = utils::inumeric(FLERR, arg[0], false, lmp);
  const char *meshfile;
  int iarg;
  if (strcmp(arg[1], "box") == 0) {
    style = BOX;
    iarg = 2;
    region = nullptr;
  } else if (strcmp(arg[1], "region") == 0) {
    style = REGION;
    if (narg < 3) utils::missing_cmd_args(FLERR, "create_atoms region", error);
    region = domain->get_region_by_id(arg[2]);
    if (!region) error->all(FLERR, "Create_atoms region {} does not exist", arg[2]);
    region->init();
    region->prematch();
    iarg = 3;
  } else if (strcmp(arg[1], "single") == 0) {
    style = SINGLE;
    if (narg < 5) utils::missing_cmd_args(FLERR, "create_atoms single", error);
    xone[0] = utils::numeric(FLERR, arg[2], false, lmp);
    xone[1] = utils::numeric(FLERR, arg[3], false, lmp);
    xone[2] = utils::numeric(FLERR, arg[4], false, lmp);
    iarg = 5;
  } else if (strcmp(arg[1], "random") == 0) {
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
  } else if (strcmp(arg[1], "mesh") == 0) {
    style = MESH;
    if (narg < 3) utils::missing_cmd_args(FLERR, "create_atoms mesh", error);
    meshfile = arg[2];
    iarg = 3;
  } else
    error->all(FLERR, "Unknown create_atoms command option {}", arg[1]);
  int scaleflag = 1;
  remapflag = 0;
  mode = ATOM;
  int molseed;
  ranmol = nullptr;
  varflag = 0;
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
  while (iarg < narg) {
    if (strcmp(arg[iarg], "basis") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms basis", error);
      int ibasis = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      int itype = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      if (ibasis <= 0 || ibasis > nbasis || itype <= 0 || itype > atom->ntypes)
        error->all(FLERR, "Out of range basis setting '{} {}' in create_atoms command", ibasis,
                   itype);
      basistype[ibasis - 1] = itype;
      iarg += 3;
    } else if (strcmp(arg[iarg], "remap") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_atoms remap", error);
      remapflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_atoms units", error);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Unknown create_atoms units option {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "var") == 0) {
      if (style == SINGLE)
        error->all(FLERR, "Cannot use 'var' keyword with 'single' style for create_atoms");
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_atoms var", error);
      delete[] vstr;
      vstr = utils::strdup(arg[iarg + 1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms set", error);
      if (strcmp(arg[iarg + 1], "x") == 0) {
        delete[] xstr;
        xstr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "y") == 0) {
        delete[] ystr;
        ystr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "z") == 0) {
        delete[] zstr;
        zstr = utils::strdup(arg[iarg + 2]);
      } else
        error->all(FLERR, "Unknown create_atoms set option {}", arg[iarg + 2]);
      iarg += 3;
    } else if (strcmp(arg[iarg], "rotate") == 0) {
      if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "create_atoms rotate", error);
      quat_user = 1;
      double thetaone;
      double axisone[3];
      thetaone = utils::numeric(FLERR, arg[iarg + 1], false, lmp) / 180.0 * MY_PI;
      axisone[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      axisone[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      axisone[2] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (axisone[0] == 0.0 && axisone[1] == 0.0 && axisone[2] == 0.0)
        error->all(FLERR, "Illegal create_atoms rotate arguments");
      if (domain->dimension == 2 && (axisone[0] != 0.0 || axisone[1] != 0.0))
        error->all(FLERR, "Invalid create_atoms rotation vector for 2d model");
      MathExtra::norm3(axisone);
      MathExtra::axisangle_to_quat(axisone, thetaone, quatone);
      iarg += 5;
    } else if (strcmp(arg[iarg], "ratio") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms ratio", error);
      subsetflag = RATIO;
      subsetfrac = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      subsetseed = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      if (subsetfrac <= 0.0 || subsetfrac > 1.0 || subsetseed <= 0)
        error->all(FLERR, "Illegal create_atoms ratio settings");
      iarg += 3;
    } else if (strcmp(arg[iarg], "subset") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms subset", error);
      subsetflag = SUBSET;
      nsubset = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      subsetseed = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      if ((nsubset <= 0) || (subsetseed <= 0))
        error->all(FLERR, "Illegal create_atoms subset settings");
      iarg += 3;
    } else if (strcmp(arg[iarg], "overlap") == 0) {
      if (style != RANDOM)
        error->all(FLERR, "Create_atoms overlap can only be used with random style");
      if (iarg + 2 > narg) error->all(FLERR, "Illegal create_atoms command");
      overlap = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (overlap <= 0) error->all(FLERR, "Illegal create_atoms overlap value: {}", overlap);
      overlapflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxtry") == 0) {
      if (style != RANDOM)
        error->all(FLERR, "Create_atoms maxtry can only be used with random style");
      if (iarg + 2 > narg) error->all(FLERR, "Illegal create_atoms command");
      maxtry = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (maxtry <= 0) error->all(FLERR, "Illegal create_atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "meshmode") == 0) {
      if (style != MESH)
        error->all(FLERR, "Create_atoms meshmode can only be used with mesh style");
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms meshmode", error);
      if (strcmp(arg[iarg + 1], "bisect") == 0) {
        mesh_style = BISECTION;
        radthresh = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      } else if (strcmp(arg[iarg + 1], "qrand") == 0) {
        mesh_style = QUASIRANDOM;
        mesh_density = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      } else
        error->all(FLERR, "Unknown create_atoms meshmode {}", arg[iarg + 2]);
      iarg += 3;
    } else if (strcmp(arg[iarg], "radscale") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_atoms radscale", error);
      if (style != MESH)
        error->all(FLERR, "Create_atoms radscale can only be used with mesh style");
      if (!atom->radius_flag)
        error->all(FLERR, "Must have atom attribute radius to set radscale factor");
      radscale = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (radscale <= 0.0) error->all(FLERR, "Illegal create_atoms radscale value: {}", radscale);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal create_atoms command option {}", arg[iarg]);
  }
  if (mode == ATOM) {
    if ((ntype <= 0) || (ntype > atom->ntypes))
      error->all(FLERR, "Invalid atom type in create_atoms command");
  }
  if (style == MESH) {
    if (mode == MOLECULE)
      error->all(FLERR, "Create_atoms mesh is not compatible with the 'mol' option");
    if (scaleflag) error->all(FLERR, "Create_atoms mesh must use 'units box' option");
  }
  ranlatt = nullptr;
  if (subsetflag != NONE) ranlatt = new RanMars(lmp, subsetseed + comm->me);
  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR, "Incomplete use of variables in create_atoms command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR, "Incomplete use of variables in create_atoms command");
  if (varflag) {
    vvar = input->variable->find(vstr);
    if (vvar < 0) error->all(FLERR, "Variable {} for create_atoms does not exist", vstr);
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR, "Variable for create_atoms is invalid style");
    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0) error->all(FLERR, "Variable {} for create_atoms does not exist", xstr);
      if (!input->variable->internalstyle(xvar))
        error->all(FLERR, "Variable for create_atoms is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0) error->all(FLERR, "Variable {} for create_atoms does not exist", ystr);
      if (!input->variable->internalstyle(yvar))
        error->all(FLERR, "Variable for create_atoms is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0) error->all(FLERR, "Variable {} for create_atoms does not exist", zstr);
      if (!input->variable->internalstyle(zvar))
        error->all(FLERR, "Variable for create_atoms is invalid style");
    }
  }
  if ((style == BOX) || (style == REGION) || (style == MESH)) {
    if (nbasis == 0) error->all(FLERR, "Cannot create atoms with undefined lattice");
  } else if (scaleflag == 1) {
    xone[0] *= domain->lattice->xlattice;
    xone[1] *= domain->lattice->ylattice;
    xone[2] *= domain->lattice->zlattice;
    overlap *= domain->lattice->xlattice;
  }
  triclinic = domain->triclinic;
  double epsilon[3];
  if (triclinic)
    epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0];
    subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1];
    subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2];
    subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0];
    subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1];
    subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2];
    subhi[2] = domain->subhi_lamda[2];
  }
  if (style == BOX || style == REGION) {
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (domain->xperiodic) {
        if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
        if (comm->myloc[0] == comm->procgrid[0] - 1) subhi[0] -= 2.0 * epsilon[0];
      }
      if (domain->yperiodic) {
        if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
        if (comm->myloc[1] == comm->procgrid[1] - 1) subhi[1] -= 2.0 * epsilon[1];
      }
      if (domain->zperiodic) {
        if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
        if (comm->myloc[2] == comm->procgrid[2] - 1) subhi[2] -= 2.0 * epsilon[2];
      }
    } else {
      if (domain->xperiodic) {
        if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
        if (comm->mysplit[0][1] == 1.0) subhi[0] -= 2.0 * epsilon[0];
      }
      if (domain->yperiodic) {
        if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
        if (comm->mysplit[1][1] == 1.0) subhi[1] -= 2.0 * epsilon[1];
      }
      if (domain->zperiodic) {
        if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
        if (comm->mysplit[2][1] == 1.0) subhi[2] -= 2.0 * epsilon[2];
      }
    }
  }
  MPI_Barrier(world);
  double time1 = platform::walltime();
  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;
  if (style == SINGLE)
    add_single();
  else if (style == RANDOM)
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
void CreateAtoms::add_single()
{
  if (remapflag) {
    imageint imagetmp = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(xone, imagetmp);
  }
  double lamda[3], *coord;
  if (triclinic) {
    domain->x2lamda(xone, lamda);
    if (remapflag) {
      if (domain->xperiodic && (lamda[0] < 0.0 || lamda[0] >= 1.0)) lamda[0] = 0.0;
      if (domain->yperiodic && (lamda[1] < 0.0 || lamda[1] >= 1.0)) lamda[1] = 0.0;
      if (domain->zperiodic && (lamda[2] < 0.0 || lamda[2] >= 1.0)) lamda[2] = 0.0;
    }
    coord = lamda;
  } else
    coord = xone;
  if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      atom->avec->create_atom(ntype, xone);
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
  if (triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    zmid = zlo + 0.5 * (zhi - zlo);
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
    zmid = zlo + 0.5 * (zhi - zlo);
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
  }
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
      if (varflag && vartest(xone) == 0) continue;
      if (triclinic) {
        domain->x2lamda(xone, lamda);
        coord = lamda;
        if (coord[0] < boxlo[0] || coord[0] >= boxhi[0] || coord[1] < boxlo[1] ||
            coord[1] >= boxhi[1] || coord[2] < boxlo[2] || coord[2] >= boxhi[2])
          continue;
      } else {
        coord = xone;
      }
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
void CreateAtoms::loop_lattice(int action)
{
  int i, j, k, m;
  const double *const *const basis = domain->lattice->basis;
  nlatt = 0;
  for (k = klo; k <= khi; k++) {
    for (j = jlo; j <= jhi; j++) {
      for (i = ilo; i <= ihi; i++) {
        for (m = 0; m < nbasis; m++) {
          double *coord;
          double x[3], lamda[3];
          x[0] = i + basis[m][0];
          x[1] = j + basis[m][1];
          x[2] = k + basis[m][2];
          domain->lattice->lattice2box(x[0], x[1], x[2]);
          if (style == REGION)
            if (!region->match(x[0], x[1], x[2])) continue;
          if (varflag && vartest(x) == 0) continue;
          if (triclinic) {
            domain->x2lamda(x, lamda);
            coord = lamda;
          } else
            coord = x;
          if (coord[0] < sublo[0] || coord[0] >= subhi[0] || coord[1] < sublo[1] ||
              coord[1] >= subhi[1] || coord[2] < sublo[2] || coord[2] >= subhi[2])
            continue;
          if (action == INSERT) {
            if (mode == ATOM) {
              atom->avec->create_atom(basistype[m], x);
            }
          } else if (action == COUNT) {
            if (nlatt == MAXSMALLINT) nlatt_overflow = 1;
          } else if (action == INSERT_SELECTED && flag[nlatt]) {
     atom->avec->create_atom(basistype[m], x);
          }
          nlatt++;
        }
      }
    }
  }
}
int CreateAtoms::vartest(double *x)
{
  if (xstr) input->variable->internal_set(xvar, x[0]);
  if (ystr) input->variable->internal_set(yvar, x[1]);
  if (zstr) input->variable->internal_set(zvar, x[2]);
  double value = input->variable->compute_equal(vvar);
  if (value == 0.0) return 0;
  return 1;
}
