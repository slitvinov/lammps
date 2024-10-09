#ifdef FIX_CLASS
// clang-format off
FixStyle(sdf_rpforce,FixSdfRepForce);
// clang-format on
#else

#ifndef LMP_FIX_SDF_REPFORCE_H
#define LMP_FIX_SDF_REPFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSdfRepForce : public Fix {
 public:
  FixSdfRepForce(class LAMMPS *, int, char **);
  ~FixSdfRepForce() override;
  int setmask() override;
  void post_force(int) override;

 private:
  float getSdf(double x, double y, double z) const;
  void getSdfGradient(double x, double y, double z, double *grad) const;

 private:
  float *sdf;
  int64_t nx_, ny_, nz_;
  double scale;
};

}    // namespace LAMMPS_NS

#endif
#endif
