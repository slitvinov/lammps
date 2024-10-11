#ifdef COMPUTE_CLASS
ComputeStyle(pe,ComputePE);
#else
#ifndef LMP_COMPUTE_PE_H
#define LMP_COMPUTE_PE_H 
#include "compute.h"
namespace LAMMPS_NS {
class ComputePE : public Compute {
 public:
  ComputePE(class LAMMPS *, int, char **);
  void init() override {}
  double compute_scalar() override;
 private:
  int pairflag, bondflag, kspaceflag, fixflag;
};
}
#endif
#endif
