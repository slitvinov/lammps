#ifdef FIX_CLASS
FixStyle(rbc, FixRBC);
#else

#ifndef LMP_FIX_RBC_H
#define LMP_FIX_RBC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRBC : public Fix {
 public:
  FixRBC(class LAMMPS *, int, char **);
  ~FixRBC() override;
  int setmask() override;
  void post_force(int) override;
  void init_list(int, class NeighList *) override;
  void init() override;
  struct RBC;

 private:
  struct RBC *rbc;
};

}    // namespace LAMMPS_NS

#endif
#endif
