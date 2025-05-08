#include <map>
#include <unordered_set>
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "region.h"
#include "region_block.h"
#include "domain.h"
using namespace LAMMPS_NS;
RegBlock::RegBlock(LAMMPS *lmp, int narg, char **arg)
    : Region(lmp, narg, arg) {
  options(narg - 8, &arg[8]);
  xlo = utils::numeric(FLERR, arg[2], false, lmp);
  xhi = utils::numeric(FLERR, arg[3], false, lmp);
  ylo = utils::numeric(FLERR, arg[4], false, lmp);
  yhi = utils::numeric(FLERR, arg[5], false, lmp);
  zlo = utils::numeric(FLERR, arg[6], false, lmp);
  zhi = utils::numeric(FLERR, arg[7], false, lmp);
}
