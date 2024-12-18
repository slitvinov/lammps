#include "region.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
Region::Region(LAMMPS *lmp, int, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr), reglist(nullptr),
      xstr(nullptr), ystr(nullptr), zstr(nullptr), tstr(nullptr) {
  id = utils::strdup(arg[0]);
  style = utils::strdup(arg[1]);
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
void Region::prematch() {}
int Region::match(double x, double y, double z) {
  return !(inside(x, y, z) ^ interior);
}
void Region::options(int narg, char **arg) {
  if (narg < 0)
    utils::missing_cmd_args(FLERR, "region", error);
  interior = 1;
  for (int i = 0; i < 6; i++)
    open_faces[i] = 0;
  xscale = yscale = zscale = 1.0;
}
void Region::reset_vel() {
  for (int i = 0; i < size_restart; i++)
    prev[i] = 0;
}
