#ifdef COMPUTE_CLASS
ComputeStyle(pressure,ComputePressure);
#else
#ifndef LMP_COMPUTE_PRESSURE_H
#define LMP_COMPUTE_PRESSURE_H 
#include "compute.h"
namespace LAMMPS_NS {
class ComputePressure : public Compute {
 public:
  ComputePressure(class LAMMPS *, int, char **);
  ~ComputePressure() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  void reset_extra_compute_fix(const char *) override;
 protected:
  double boltz, nktv2p, inv_volume;
  int nvirial, dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag, pairflag, bondflag;
  int fixflag, kspaceflag;
  void virial_compute(int, int);
 private:
  char *pstyle;
  int nsub;
};
}
#endif
#endif
