#ifdef FIX_CLASS
FixStyle(READ_RESTART,FixReadRestart);
#else
#ifndef LMP_FIX_READ_RESTART_H
#define LMP_FIX_READ_RESTART_H 
#include "fix.h"
namespace LAMMPS_NS {
class FixReadRestart : public Fix {
 public:
  int *count;
  double **extra;
  FixReadRestart(class LAMMPS *, int, char **);
  ~FixReadRestart() override;
  int setmask() override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
 private:
  int nextra;
};
}
#endif
#endif
