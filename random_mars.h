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
  double gaussian(double mu, double sigma);
  double rayleigh(double sigma);
  double besselexp(double theta, double alpha, double cp);
  void select_subset(bigint, int, int *, int *);
  void get_state(double *);
  void set_state(double *);

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
