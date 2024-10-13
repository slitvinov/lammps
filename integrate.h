#ifndef LMP_INTEGRATE_H
#define LMP_INTEGRATE_H 
#include "pointers.h"
namespace LAMMPS_NS {
class Integrate : protected Pointers {
 public:
  Integrate(class LAMMPS *, int, char **);
  ~Integrate() override;
  virtual void init();
  virtual void setup(int flag) = 0;
  virtual void run(int) = 0;
  virtual void force_clear() = 0;
  virtual void cleanup() {}
  virtual void reset_dt() {}
 protected:
  int eflag, vflag;
  int virial_style;
  int external_force_clear;
  int nelist_global, nelist_atom;
  int nvlist_global, nvlist_atom, ncvlist_atom;
  class Compute **elist_global;
  class Compute **elist_atom;
  class Compute **vlist_global;
  class Compute **vlist_atom;
  class Compute **cvlist_atom;
  int pair_compute_flag;
  int kspace_compute_flag;
  void ev_setup();
  void ev_set(bigint);
};
}
#endif
