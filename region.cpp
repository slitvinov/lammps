#include <unordered_set>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "region.h"
#include "domain.h"
#include "lmptype.h"
#include "update.h"
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
void Region::init() { vel_timestep = -1; }
void Region::prematch() {}
int Region::match(double x, double y, double z) {
  return !(inside(x, y, z) ^ interior);
}
void Region::options(int narg, char **arg) {
  interior = 1;
  for (int i = 0; i < 6; i++)
    open_faces[i] = 0;
  xscale = yscale = zscale = 1.0;
}
void Region::reset_vel() {
  for (int i = 0; i < size_restart; i++)
    prev[i] = 0;
}
