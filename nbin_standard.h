#ifdef NBIN_CLASS
NBinStyle(standard,
          NBinStandard,
          NB_STANDARD);
#else
#ifndef LMP_NBIN_STANDARD_H
#define LMP_NBIN_STANDARD_H 
#include "nbin.h"
namespace LAMMPS_NS {
class NBinStandard : public NBin {
 public:
  NBinStandard(class LAMMPS *);
  void bin_atoms_setup(int) override;
  void setup_bins(int) override;
  void bin_atoms() override;
  double memory_usage() override;
};
}
#endif
#endif
