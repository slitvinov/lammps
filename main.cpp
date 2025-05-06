#include <set>
#include <map>
#include <unordered_set>
#include "pointers.h"
#include "atom.h"
#include "atom_vec_atomic.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "error.h"
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
#include "create_atoms.h"
#include "run.h"
#include "universe.h"
#include "update.h"
#include "verlet.h"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#define LAMMPS_VERSION "28 Mar 2023"
#define LAMMPS_UPDATE "Development"
#if defined(LAMMPS_UPDATE)
#define UPDATE_STRING " - " LAMMPS_UPDATE
#else
#define UPDATE_STRING ""
#endif
static void print_style(FILE *fp, const char *str, int &pos);
using namespace LAMMPS_NS;
LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator)
    : memory(nullptr), error(nullptr), universe(nullptr), input(nullptr),
      atom(nullptr), update(nullptr), neighbor(nullptr), comm(nullptr),
      domain(nullptr), force(nullptr), modify(nullptr), group(nullptr) {
  memory = new Memory(this);
  error = new Error(this);
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
    if (strcmp(arg[iarg], "-in") == 0 || strcmp(arg[iarg], "-i") == 0) {
      inflag = iarg + 1;
      iarg += 2;
    }
  }
  if (universe->existflag == 0)
    universe->add_world(nullptr);
  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag], "none") == 0)
      universe->uscreen = nullptr;
    else {
      universe->uscreen = fopen(arg[screenflag], "w");
      if (universe->uscreen == nullptr)
        error->universe_one(
            FLERR, fmt::format("Cannot open universe screen file {}: {}",
                               arg[screenflag], utils::getsyserror()));
    }
    if (logflag == 0) {
      if (helpflag == 0) {
        universe->ulogfile = fopen("log.lammps", "w");
        if (universe->ulogfile == nullptr)
          error->universe_warn(FLERR, "Cannot open log.lammps for writing: " +
                                          utils::getsyserror());
      }
    } else if (strcmp(arg[logflag], "none") == 0)
      universe->ulogfile = nullptr;
    else {
      universe->ulogfile = fopen(arg[logflag], "w");
      if (universe->ulogfile == nullptr)
        error->universe_one(FLERR,
                            fmt::format("Cannot open universe log file {}: {}",
                                        arg[logflag], utils::getsyserror()));
    }
  }
  if (universe->me > 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else
      universe->uscreen = nullptr;
    universe->ulogfile = nullptr;
  }
  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;
    if (universe->me == 0) {
      if (inflag == 0)
        infile = stdin;
      else if (strcmp(arg[inflag], "none") == 0)
        infile = stdin;
      else
        infile = fopen(arg[inflag], "r");
      if (infile == nullptr)
        error->one(FLERR, "Cannot open input script {}: {}", arg[inflag],
                   utils::getsyserror());
      if ((inflag == 0) && (universe->nprocs > 1))
        error->warning(
            FLERR, "Using I/O redirection is unreliable with parallel runs. "
                   "Better to use the -in switch to read input files.");
      utils::flush_buffers(this);
    }
  }
  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR, "Smallint setting in lmptype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR, "Imageint setting in lmptype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR, "Tagint setting in lmptype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR, "Bigint setting in lmptype.h is invalid");
  int mpisize;
  MPI_Type_size(MPI_LMP_TAGINT, &mpisize);
  if (mpisize != sizeof(tagint))
    error->all(FLERR,
               "MPI_LMP_TAGINT and tagint in lmptype.h are not compatible");
  MPI_Type_size(MPI_LMP_BIGINT, &mpisize);
  if (mpisize != sizeof(bigint))
    error->all(FLERR,
               "MPI_LMP_BIGINT and bigint in lmptype.h are not compatible");
#ifdef LAMMPS_SMALLBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || sizeof(tagint) != 4 ||
      sizeof(bigint) != 8)
    error->all(FLERR, "Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_BIGBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 8 || sizeof(tagint) != 8 ||
      sizeof(bigint) != 8)
    error->all(FLERR, "Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_SMALLSMALL
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || sizeof(tagint) != 4 ||
      sizeof(bigint) != 4)
    error->all(FLERR, "Small to big integers are not sized correctly");
#endif
  input = new Input(this, narg, arg);
  if (helpflag) {
    error->done(0);
  } else {
    create();
  }
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
