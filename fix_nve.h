#ifdef FIX_CLASS
FixStyle(nve, FixNVE);
#else
#ifndef LMP_FIX_NVE_H
#define LMP_FIX_NVE_H
#include "fix.h"
namespace LAMMPS_NS {
class FixNVE : public Fix {
public:
  FixNVE(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;
};
} // namespace LAMMPS_NS
#endif
#endif
