#ifdef FIX_CLASS
// clang-format off
FixStyle(fileforce,FixFileforce);
// clang-format on
#else

#ifndef LMP_FIX_FILEFORCE_H
#define LMP_FIX_FILEFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFileforce : public Fix {
 public:
  FixFileforce(class LAMMPS *, int, char **);
  ~FixFileforce() override;
  int setmask() override;
  void post_force(int) override;
  struct Force;

 private:
  struct Force *force;
};

}    // namespace LAMMPS_NS

#endif
#endif
