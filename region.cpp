#include "region.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
Region::Region(LAMMPS *lmp, int, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr), reglist(nullptr),
      contact(nullptr), xstr(nullptr), ystr(nullptr), zstr(nullptr),
      tstr(nullptr) {
  id = utils::strdup(arg[0]);
  style = utils::strdup(arg[1]);
  varshape = 0;
  xstr = ystr = zstr = tstr = nullptr;
  dx = dy = dz = 0.0;
  size_restart = 5;
  Region::reset_vel();
  nregion = 1;
}
Region::~Region() {
  delete[] id;
  delete[] style;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] tstr;
}
void Region::init() { vel_timestep = -1; }
void Region::prematch() {
  if (varshape)
    shape_update();
}
int Region::match(double x, double y, double z) {
  if (openflag)
    return 1;
  return !(inside(x, y, z) ^ interior);
}
void Region::options(int narg, char **arg) {
  if (narg < 0)
    utils::missing_cmd_args(FLERR, "region", error);
  interior = 1;
  scaleflag = 1;
  moveflag = 0;
  openflag = 0;
  for (int i = 0; i < 6; i++)
    open_faces[i] = 0;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "region units", error);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Illegal region units: {}", arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal region command argument: {}", arg[iarg]);
  }
  if ((moveflag) &&
      (strcmp(style, "union") == 0 || strcmp(style, "intersect") == 0))
    error->all(FLERR, "Region union or intersect cannot be dynamic");
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else
    xscale = yscale = zscale = 1.0;
  if (moveflag)
    dynamic = 1;
  else
    dynamic = 0;
}
void Region::reset_vel() {
  for (int i = 0; i < size_restart; i++)
    prev[i] = 0;
}
