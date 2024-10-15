#ifndef LMP_RANMARS_H
#define LMP_RANMARS_H
#include "pointers.h"
namespace LAMMPS_NS {
class RanMars : protected Pointers {
public:
  RanMars(class LAMMPS *, int);
  ~RanMars() override;
  double uniform();
  double gaussian();
private:
  char padding[1024];
  int save;
  double second;
  double *u;
  int i97, j97;
  double c, cd, cm;
};
} // namespace LAMMPS_NS
#endif
