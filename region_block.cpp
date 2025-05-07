#include <map>
#include <unordered_set>
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "region.h"
#include "region_block.h"
#include "domain.h"
#include "input.h"
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
  face[0][0] = -1.0;
  face[0][1] = 0.0;
  face[0][2] = 0.0;
  face[1][0] = 1.0;
  face[1][1] = 0.0;
  face[1][2] = 0.0;
  face[2][0] = 0.0;
  face[2][1] = -1.0;
  face[2][2] = 0.0;
  face[3][0] = 0.0;
  face[3][1] = 1.0;
  face[3][2] = 0.0;
  face[4][0] = 0.0;
  face[4][1] = 0.0;
  face[4][2] = -1.0;
  face[5][0] = 0.0;
  face[5][1] = 0.0;
  face[5][2] = 1.0;
  corners[0][0][0] = xlo;
  corners[0][0][1] = ylo;
  corners[0][0][2] = zlo;
  corners[0][1][0] = xlo;
  corners[0][1][1] = ylo;
  corners[0][1][2] = zhi;
  corners[0][2][0] = xlo;
  corners[0][2][1] = yhi;
  corners[0][2][2] = zhi;
  corners[0][3][0] = xlo;
  corners[0][3][1] = yhi;
  corners[0][3][2] = zlo;
  corners[1][0][0] = xhi;
  corners[1][0][1] = ylo;
  corners[1][0][2] = zlo;
  corners[1][1][0] = xhi;
  corners[1][1][1] = ylo;
  corners[1][1][2] = zhi;
  corners[1][2][0] = xhi;
  corners[1][2][1] = yhi;
  corners[1][2][2] = zhi;
  corners[1][3][0] = xhi;
  corners[1][3][1] = yhi;
  corners[1][3][2] = zlo;
  copy3(corners[0][0], corners[2][0]);
  copy3(corners[1][0], corners[2][1]);
  copy3(corners[1][1], corners[2][2]);
  copy3(corners[0][1], corners[2][3]);
  copy3(corners[0][3], corners[3][0]);
  copy3(corners[0][2], corners[3][1]);
  copy3(corners[1][2], corners[3][2]);
  copy3(corners[1][3], corners[3][3]);
  copy3(corners[0][0], corners[4][0]);
  copy3(corners[0][3], corners[4][1]);
  copy3(corners[1][3], corners[4][2]);
  copy3(corners[1][0], corners[4][3]);
  copy3(corners[0][1], corners[5][0]);
  copy3(corners[1][1], corners[5][1]);
  copy3(corners[1][2], corners[5][2]);
  copy3(corners[0][2], corners[5][3]);
}
void RegBlock::init() { Region::init(); }
int RegBlock::inside(double x, double y, double z) {
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}
