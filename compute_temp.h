#ifdef COMPUTE_CLASS
ComputeStyle(temp,ComputeTemp);
#else
#ifndef LMP_COMPUTE_TEMP_H
#define LMP_COMPUTE_TEMP_H 
#include "compute.h"
namespace LAMMPS_NS {
class ComputeTemp : public Compute {
 public:
  ComputeTemp(class LAMMPS *, int, char **);
  ~ComputeTemp() override;
  void init() override {}
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;
 protected:
  double tfactor;
  virtual void dof_compute();
};
}
#endif
#endif
