#ifdef FIX_CLASS
FixStyle(pair,FixPair);
#else
#ifndef LMP_FIX_PAIR_H
#define LMP_FIX_PAIR_H 
#include "fix.h"
namespace LAMMPS_NS {
class FixPair : public Fix {
 public:
  FixPair(class LAMMPS *, int, char **);
  ~FixPair() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void min_pre_force(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  double memory_usage() override;
 private:
  int nevery, nfield, ncols;
  bigint lasttime;
  char *pairname;
  char **fieldname, **triggername;
  int *trigger;
  int **triggerptr;
  class Pair *pstyle;
  double *vector;
  double **array;
  void query_pstyle(LAMMPS *lmp);
};
}
#endif
#endif
