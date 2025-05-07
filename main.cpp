#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <map>
#include <mpi.h>
#include <set>
#include <stdlib.h>
#include <unistd.h>
#include <unordered_set>
#include <string>
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "atom.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair_dpd.h"
#include "region_block.h"
#include "universe.h"
#include "update.h"
using namespace LAMMPS_NS;
LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator)
    : memory(nullptr), error(nullptr), universe(nullptr), input(nullptr),
      atom(nullptr), update(nullptr), neighbor(nullptr), comm(nullptr),
      domain(nullptr), force(nullptr), modify(nullptr), group(nullptr) {
  memory = new Memory(this);
  universe = new Universe(this, communicator);
  restart_ver = -1;
  external_comm = 0;
  skiprunflag = 0;
  screen = nullptr;
  logfile = nullptr;
  infile = nullptr;
  int iarg = 1;
  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int helpflag = 0;
  int nonbufflag = 0;
  pair_only_flag = 0;
  exename = arg[0];
  packargs = nullptr;
  num_package = 0;
  char *restartfile = nullptr;
  int wfirst, wlast;
  int kkfirst, kklast;
  int npack = 0;
  int *pfirst = nullptr;
  int *plast = nullptr;
  iarg = 1;
  while (iarg < narg) {
    inflag = iarg + 1;
    iarg += 2;
  }
  universe->add_world(nullptr);
  if (universe->me == 0) {
    universe->uscreen = stdout;
    universe->ulogfile = fopen("log.lammps", "w");
  }
  if (universe->me > 0) {
    universe->uscreen = stdout;
    universe->ulogfile = nullptr;
  }
  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;
    if (universe->me == 0) {
      infile = fopen(arg[inflag], "r");
    }
  }
  int mpisize;
  MPI_Type_size(MPI_LMP_TAGINT, &mpisize);
  MPI_Type_size(MPI_LMP_BIGINT, &mpisize);
  input = new Input(this, narg, arg);
  create();
}
void LAMMPS::create() {
  force = nullptr;
  comm = new CommBrick(this);
  neighbor = new Neighbor(this);
  domain = new Domain(this);
  atom = new Atom(this);
  atom->create_avec("atomic", 0, nullptr, 1);
  group = new Group(this);
  force = new Force(this);
  modify = new Modify(this);
  update = new Update(this);
}
void LAMMPS::init() {
  update->init();
  force->init();
  domain->init();
  atom->init();
  modify->init();
  neighbor->init();
  comm->init();
}
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm lammps_comm = MPI_COMM_WORLD;
  auto lammps = new LAMMPS(argc, argv, lammps_comm);
  lammps->input->file();
  delete lammps;
  MPI_Barrier(lammps_comm);
  MPI_Finalize();
}
