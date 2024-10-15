#include "region_block.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_extra.h"
#include <cstring>
using namespace LAMMPS_NS;
enum { CONSTANT, VARIABLE };
#define BIG 1.0e20
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
  if (xlo > xhi)
    error->all(FLERR, "Illegal region block xlo: {} >= xhi: {}", xlo, xhi);
  if (ylo > yhi)
    error->all(FLERR, "Illegal region block ylo: {} >= yhi: {}", ylo, yhi);
  if (zlo > zhi)
    error->all(FLERR, "Illegal region block zlo: {} >= zhi: {}", zlo, zhi);
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
  contact = new Contact[cmax];
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
  MathExtra::copy3(corners[0][0], corners[2][0]);
  MathExtra::copy3(corners[1][0], corners[2][1]);
  MathExtra::copy3(corners[1][1], corners[2][2]);
  MathExtra::copy3(corners[0][1], corners[2][3]);
  MathExtra::copy3(corners[0][3], corners[3][0]);
  MathExtra::copy3(corners[0][2], corners[3][1]);
  MathExtra::copy3(corners[1][2], corners[3][2]);
  MathExtra::copy3(corners[1][3], corners[3][3]);
  MathExtra::copy3(corners[0][0], corners[4][0]);
  MathExtra::copy3(corners[0][3], corners[4][1]);
  MathExtra::copy3(corners[1][3], corners[4][2]);
  MathExtra::copy3(corners[1][0], corners[4][3]);
  MathExtra::copy3(corners[0][1], corners[5][0]);
  MathExtra::copy3(corners[1][1], corners[5][1]);
  MathExtra::copy3(corners[1][2], corners[5][2]);
  MathExtra::copy3(corners[0][2], corners[5][3]);
}
RegBlock::~RegBlock() {
  delete[] xlostr;
  delete[] xhistr;
  delete[] ylostr;
  delete[] yhistr;
  delete[] zlostr;
  delete[] zhistr;
  delete[] contact;
}
void RegBlock::init() { Region::init(); }
int RegBlock::inside(double x, double y, double z) {
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}
