#ifdef FIX_CLASS
FixStyle(fileforce,FixFileforce);
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
}
#endif
#endif
