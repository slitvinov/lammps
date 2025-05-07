#include "universe.h"
#include "memory.h"
#include <cstring>
using namespace LAMMPS_NS;
#define MAXLINE 256
Universe::Universe(LAMMPS *lmp, MPI_Comm communicator) : Pointers(lmp) {
  uworld = uorig = communicator;
  MPI_Comm_rank(uworld, &me);
  MPI_Comm_size(uworld, &nprocs);
  uscreen = stdout;
  ulogfile = nullptr;
  existflag = 0;
  nworlds = 0;
  procs_per_world = nullptr;
  root_proc = nullptr;
  memory->create(uni2orig, nprocs, "universe:uni2orig");
  for (int i = 0; i < nprocs; i++)
    uni2orig[i] = i;
}
void Universe::add_world(char *str) {
  int n, nper;
  n = 1;
  nper = 0;
  nper = nprocs;
  memory->grow(procs_per_world, nworlds + n, "universe:procs_per_world");
  memory->grow(root_proc, (nworlds + n), "universe:root_proc");
  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    root_proc[nworlds] = 0;
    if (me >= root_proc[nworlds])
      iworld = nworlds;
    nworlds++;
  }
}
