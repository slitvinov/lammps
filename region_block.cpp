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
enum { CONSTANT, VARIABLE };
#define BIG 1.0e20
static void copy3(const double *v, double *ans) {
  ans[0] = v[0];
  ans[1] = v[1];
  ans[2] = v[2];
}

RegBlock::RegBlock(LAMMPS *lmp, int narg, char **arg)
    : Region(lmp, narg, arg), xlostr(nullptr), xhistr(nullptr), ylostr(nullptr),
      yhistr(nullptr), zlostr(nullptr), zhistr(nullptr) {
  options(narg - 8, &arg[8]);
  xlostyle = CONSTANT;
  xlo = xscale * utils::numeric(FLERR, arg[2], false, lmp);
  xhistyle = CONSTANT;
  xhi = xscale * utils::numeric(FLERR, arg[3], false, lmp);
  ylostyle = CONSTANT;
  ylo = yscale * utils::numeric(FLERR, arg[4], false, lmp);
  yhistyle = CONSTANT;
  yhi = yscale * utils::numeric(FLERR, arg[5], false, lmp);
  zlostyle = CONSTANT;
  zlo = zscale * utils::numeric(FLERR, arg[6], false, lmp);
  zhistyle = CONSTANT;
  zhi = zscale * utils::numeric(FLERR, arg[7], false, lmp);
  if (interior) {
    bboxflag = 1;
    extent_xlo = xlo;
    extent_xhi = xhi;
    extent_ylo = ylo;
    extent_yhi = yhi;
    extent_zlo = zlo;
    extent_zhi = zhi;
  } else
    bboxflag = 0;
  cmax = 6;
  if (interior)
    tmax = 3;
  else
    tmax = 1;
}
void RegBlock::init() { Region::init(); }
int RegBlock::inside(double x, double y, double z) {
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}
