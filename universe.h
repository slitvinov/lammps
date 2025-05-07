#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H
namespace LAMMPS_NS {
class Universe : protected Pointers {
public:
  MPI_Comm uworld;
  int me, nprocs;
  FILE *uscreen;
  FILE *ulogfile;
  int existflag;
  int nworlds;
  int iworld;
  int *procs_per_world;
  int *root_proc;
  MPI_Comm uorig;
  int *uni2orig;
  Universe(class LAMMPS *, MPI_Comm);
  void add_world(char *);
};
} // namespace LAMMPS_NS
#endif
