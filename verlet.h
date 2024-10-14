#ifdef INTEGRATE_CLASS
IntegrateStyle(verlet, Verlet);
#else
#ifndef LMP_VERLET_H
#define LMP_VERLET_H
#include "integrate.h"
namespace LAMMPS_NS {
class Verlet : public Integrate {
public:
  Verlet(class LAMMPS *, int, char **);
  void init() override;
  void setup(int flag) override;
  void run(int) override;
  void force_clear() override;
  void cleanup() override;

protected:
  int triclinic;
  int torqueflag, extraflag;
};
} // namespace LAMMPS_NS
#endif
#endif
