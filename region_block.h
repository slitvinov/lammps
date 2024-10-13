#ifdef REGION_CLASS
RegionStyle(block,RegBlock);
#else
#ifndef LMP_REGION_BLOCK_H
#define LMP_REGION_BLOCK_H 
#include "region.h"
namespace LAMMPS_NS {
class RegBlock : public Region {
  friend class FixPour;
 public:
  RegBlock(class LAMMPS *, int, char **);
  ~RegBlock() override;
  void init() override;
  int inside(double, double, double) override;
 protected:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double corners[6][4][3];
  double face[6][3];
  int xlostyle, xlovar, xhistyle, xhivar;
  int ylostyle, ylovar, yhistyle, yhivar;
  int zlostyle, zlovar, zhistyle, zhivar;
  char *xlostr, *ylostr, *zlostr;
  char *xhistr, *yhistr, *zhistr;
};
}
#endif
#endif
