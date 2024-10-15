#ifndef LMP_RANPARK_H
#define LMP_RANPARK_H
#include "pointers.h"
namespace LAMMPS_NS {
class RanPark : protected Pointers {
public:
  RanPark(class LAMMPS *, int);
  double uniform();

private:
  int seed, save;
  double second;
};
} // namespace LAMMPS_NS
#endif
