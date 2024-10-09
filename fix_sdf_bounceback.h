#ifdef FIX_CLASS
// clang-format off
FixStyle(sdf_bounceback,FixSdfBounceback);
// clang-format on
#else

#ifndef LMP_FIX_SDF_BOUNCEBACK_H
#define LMP_FIX_SDF_BOUNCEBACK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSdfBounceback : public Fix {
 public:
  FixSdfBounceback(class LAMMPS *, int, char **);
  ~FixSdfBounceback() override;
  int setmask() override;
  void post_integrate() override;

 private:
  float getSdf(double x, double y, double z) const;

 private:
  float *sdf;
  int64_t nx_, ny_, nz_;
};

}    // namespace LAMMPS_NS

#endif
#endif
