#include <cmath>
#include <cstring>
#include <map>
#include <unordered_set>
#include <string>
#include <vector>
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "lattice.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "update.h"
using namespace LAMMPS_NS;
Lattice::Lattice(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp) {
  nbasis = 0;
  basis = nullptr;
  style = NONE;
  if (style == NONE) {
    xlattice = ylattice = zlattice = utils::numeric(FLERR, arg[1], false, lmp);
    return;
  }
}
