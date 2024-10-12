#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H 
#include <cstdio>
#include <mpi.h>
namespace LAMMPS_NS {
class LAMMPS {
 public:
  class Memory *memory;
  class Error *error;
  class Universe *universe;
  class Input *input;
  class Atom *atom;
  class Update *update;
  class Neighbor *neighbor;
  class Comm *comm;
  class Domain *domain;
  class Force *force;
  class Modify *modify;
  class Group *group;
  int restart_ver;
  MPI_Comm world;
  FILE *infile;
  FILE *screen;
  FILE *logfile;
  double initclock;
  int skiprunflag;
  int pair_only_flag;
  char *exename;
  char ***packargs;
  int num_package;
  MPI_Comm external_comm;
  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void init();
  void destroy();
 private:
  LAMMPS(){};
  LAMMPS(const LAMMPS &){};
};
}
#endif
