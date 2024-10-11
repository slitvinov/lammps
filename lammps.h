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
  const char *version;
  int num_ver;
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
  static bool has_git_info();
  static const char *git_commit();
  static const char *git_branch();
  static const char *git_descriptor();
  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void post_create();
  void init();
  void destroy();
 private:
  LAMMPS(){};
  LAMMPS(const LAMMPS &){};
};
}
#endif
