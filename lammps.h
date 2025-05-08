#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H
namespace LAMMPS_NS {
class LAMMPS {
public:
  class Memory *memory;
  class Universe *universe;
  class Input *input;
  class Atom *atom;
  class Update *update;
  class Neighbor *neighbor;
  class Comm *comm;
  class Domain *domain;
  class Force *force;
  class Group *group;
  int restart_ver;
  MPI_Comm world;
  FILE *infile;
  FILE *screen;
  FILE *logfile;
  int skiprunflag;
  int pair_only_flag;
  char *exename;
  char ***packargs;
  int num_package;
  MPI_Comm external_comm;
  LAMMPS(int, char **, MPI_Comm);
  void create();
  void init();

  class FixNVE *fix_nve;
  class RegBlock *region_block;

private:
  LAMMPS(){};
  LAMMPS(const LAMMPS &){};
};
} // namespace LAMMPS_NS
#endif
