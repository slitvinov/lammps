#ifndef LMP_UPDATE_H
#define LMP_UPDATE_H 
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Update : protected Pointers {
 public:
  double dt;
  double etol, ftol;
  bigint ntimestep;
  int nsteps;
  int whichflag;
  double atime;
  bigint atimestep;
  bigint firststep, laststep;
  bigint beginstep, endstep;
  int first_update;
  int max_eval;
  int restrict_output;
  int setupflag;
  int multireplica;
  int dt_default;
  bigint eflag_global, eflag_atom;
  bigint vflag_global, vflag_atom;
  char *unit_style;
  class Integrate *integrate;
  char *integrate_style;
  typedef Integrate *(*IntegrateCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, IntegrateCreator> IntegrateCreatorMap;
  IntegrateCreatorMap *integrate_map;
  Update(class LAMMPS *);
  ~Update() override;
  void init();
  void set_units(const char *);
  void create_integrate(int, char **, int);
  void reset_timestep(int, char **);
  void reset_timestep(bigint, bool);
  void update_time();
  double memory_usage();
 private:
  void new_integrate(char *, int, char **);
};
}
#endif
