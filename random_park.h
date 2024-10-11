#ifndef LMP_RANPARK_H
#define LMP_RANPARK_H 
#include "pointers.h"
namespace LAMMPS_NS {
class RanPark : protected Pointers {
 public:
  RanPark(class LAMMPS *, int);
  double uniform();
  double gaussian();
  void reset(int);
  void reset(int, double *);
  int state();
 private:
  int seed, save;
  double second;
};
}
#endif
