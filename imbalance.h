#ifndef LMP_IMBALANCE_H
#define LMP_IMBALANCE_H 
#include "pointers.h"
#include <string>
namespace LAMMPS_NS {
class Imbalance : protected Pointers {
 public:
  Imbalance(class LAMMPS *);
  virtual int options(int, char **) = 0;
  virtual void init(int){};
  virtual void compute(double *) = 0;
  virtual std::string info() = 0;
};
}
#endif
