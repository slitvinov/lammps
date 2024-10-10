// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"

#include "style_atom.h"      // IWYU pragma: keep
#include "style_command.h"   // IWYU pragma: keep
#include "style_integrate.h" // IWYU pragma: keep
#include "style_pair.h"      // IWYU pragma: keep
#include "style_region.h"    // IWYU pragma: keep

#include "atom.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "suffix.h"
#include "universe.h"
#include "update.h"
#include "variable.h"
#include "version.h"

#if defined(LMP_PLUGIN)
#include "plugin.h"
#endif

#include <cctype>
#include <cmath>
#include <cstring>
#include <map>
#include "lmpgitversion.h"

#if defined(LAMMPS_UPDATE)
#define UPDATE_STRING " - " LAMMPS_UPDATE
#else
#define UPDATE_STRING ""
#endif

static void print_style(FILE *fp, const char *str, int &pos);
using namespace LAMMPS_NS;

/** \class LAMMPS_NS::LAMMPS
 * \brief LAMMPS simulation instance
 *
 * The LAMMPS class contains pointers of all constituent class instances
 * and global variables that are used by a LAMMPS simulation. Its contents
 * represent the entire state of the simulation.
 *
 * The LAMMPS class manages the components of an MD simulation by creating,
 * deleting, and initializing instances of the classes it is composed of,
 * processing command line flags, and providing access to some global properties.
 * The specifics of setting up and running a simulation are handled by the
 * individual component class instances. */

/** Create a LAMMPS simulation instance
 *
 * The LAMMPS constructor starts up a simulation by allocating all
 * fundamental classes in the necessary order, parses input switches
 * and their arguments, initializes communicators, screen and logfile
 * output FILE pointers.
 *
 * \param narg number of arguments
 * \param arg list of arguments
 * \param communicator MPI communicator used by this LAMMPS instance
 */
LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator) :
  memory(nullptr), error(nullptr), universe(nullptr), input(nullptr), atom(nullptr),
  update(nullptr), neighbor(nullptr), comm(nullptr), domain(nullptr), force(nullptr),
  modify(nullptr), group(nullptr)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);

  version = (const char *) LAMMPS_VERSION;
  num_ver = utils::date2num(version);
  restart_ver = -1;

  external_comm = 0;
  skiprunflag = 0;

  screen = nullptr;
  logfile = nullptr;
  infile = nullptr;

  initclock = platform::walltime();

#if defined(LMP_PYTHON) && defined(_WIN32)
  // If the LAMMPSHOME environment variable is set, it should point
  // to the location of the LAMMPS installation tree where we bundle
  // the matching Python installation for use with the PYTHON package.
  // This is currently only used on Windows with the Windows installer packages
  const char *lmpenv = getenv("LAMMPSHOME");
  if (lmpenv) {
    platform::putenv(fmt::format("PYTHONHOME={}",lmpenv));
  }
#endif

  // check if -mpicolor is first arg
  // if so, then 2 or more apps were launched with one mpirun command
  //   this means passed communicator (e.g. MPI_COMM_WORLD) is bigger than LAMMPS
  //   universe communicator needs to shrink to be just LAMMPS
  // syntax: -mpicolor color
  //   color = integer for this app, different than any other app(s)
  // do the following:
  //   perform an MPI_Comm_split() to create a new LAMMPS-only subcomm
  //   NOTE: this assumes other app(s) make same call, else will hang!
  //   re-create universe with subcomm
  //   store comm that all apps belong to in external_comm

  int iarg = 1;
  if (narg-iarg >= 2 && (strcmp(arg[iarg],"-mpicolor") == 0 ||
                         strcmp(arg[iarg],"-m") == 0)) {
    int me,nprocs;
    MPI_Comm_rank(communicator,&me);
    MPI_Comm_size(communicator,&nprocs);
    int color = atoi(arg[iarg+1]);
    MPI_Comm subcomm;
    MPI_Comm_split(communicator,color,me,&subcomm);
    external_comm = communicator;
    communicator = subcomm;
    delete universe;
    universe = new Universe(this,communicator);
  }

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int restart2data = 0;
  int restart2dump = 0;
  int restartremap = 0;
  int helpflag = 0;
  int nonbufflag = 0;

  pair_only_flag = 0;
  if (arg) exename = arg[0];
  else exename = nullptr;
  packargs = nullptr;
  num_package = 0;
  char *restartfile = nullptr;
  int wfirst,wlast;
  int kkfirst,kklast;

  int npack = 0;
  int *pfirst = nullptr;
  int *plast = nullptr;

  iarg = 1;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;

    } else if (strcmp(arg[iarg],"-in") == 0 ||
               strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-mpicolor") == 0 ||
               strcmp(arg[iarg],"-m") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (iarg != 1) error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-nonbuf") == 0 ||
               strcmp(arg[iarg],"-nb") == 0) {
      nonbufflag = 1;
      iarg++;

    } else if (strcmp(arg[iarg],"-package") == 0 ||
               strcmp(arg[iarg],"-pk") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      memory->grow(pfirst,npack+1,"lammps:pfirst");
      memory->grow(plast,npack+1,"lammps:plast");
      // delimit args for package command invocation
      // any package arg with leading "-" will be followed by numeric digit
      iarg++;
      pfirst[npack] = iarg;
      while (iarg < narg) {
        if (arg[iarg][0] != '-') iarg++;
        else if (isdigit(arg[iarg][1])) iarg++;
        else break;
      }
      plast[npack++] = iarg;

    } else if (strcmp(arg[iarg],"-partition") == 0 ||
        strcmp(arg[iarg],"-p") == 0) {
      universe->existflag = 1;
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
        universe->add_world(arg[iarg]);
        iarg++;
      }

    } else if (strcmp(arg[iarg],"-plog") == 0 ||
               strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
               strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-reorder") == 0 ||
               strcmp(arg[iarg],"-ro") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (universe->existflag)
        error->universe_all(FLERR,"Cannot use -reorder after -partition");
      universe->reorder(arg[iarg+1],arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"-restart2data") == 0 ||
               strcmp(arg[iarg],"-r2data") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (restart2dump)
        error->universe_all(FLERR,
                            "Cannot use both -restart2data and -restart2dump");
      restart2data = 1;
      restartfile = arg[iarg+1];
      // check for restart remap flag
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
      // delimit args for the write_data command
      wfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wlast = iarg;

    } else if (strcmp(arg[iarg],"-restart2dump") == 0 ||
               strcmp(arg[iarg],"-r2dump") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (restart2data)
        error->universe_all(FLERR,
                            "Cannot use both -restart2data and -restart2dump");
      restart2dump = 1;
      restartfile = arg[iarg+1];
      // check for restart remap flag
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
      // delimit args for the write_dump command
      wfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wlast = iarg;

    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-skiprun") == 0 ||
               strcmp(arg[iarg],"-sr") == 0) {
      skiprunflag = 1;
      ++iarg;
    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;

    } else {
      error->universe_all(FLERR, fmt::format("Invalid command-line argument: {}", arg[iarg]) );
    }
  }

  // if no partition command-line switch, universe is one world with all procs

  if (universe->existflag == 0) universe->add_world(nullptr);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions do not match number of allocated processors");

  // universe cannot use stdin for input file

  if (universe->existflag && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");

  // if no partition command-line switch, cannot use -pscreen option

  if (universe->existflag == 0 && partscreenflag)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");

  // if no partition command-line switch, cannot use -plog option

  if (universe->existflag == 0 && partlogflag)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = nullptr;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open universe screen file {}: {}",
                                              arg[screenflag],utils::getsyserror()));
    }
    if (logflag == 0) {
      if (helpflag == 0) {
        universe->ulogfile = fopen("log.lammps","w");
        if (universe->ulogfile == nullptr)
          error->universe_warn(FLERR,"Cannot open log.lammps for writing: "
                               + utils::getsyserror());
      }
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = nullptr;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open universe log file {}: {}",
                                              arg[logflag],utils::getsyserror()));
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = nullptr;
    universe->ulogfile = nullptr;
  }

  // make universe and single world the same, since no partition switch
  // world inherits settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else if (strcmp(arg[inflag], "none") == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == nullptr)
        error->one(FLERR,"Cannot open input script {}: {}", arg[inflag], utils::getsyserror());
      if (!helpflag)
        utils::logmesg(this,fmt::format("LAMMPS ({}{})\n",version,UPDATE_STRING));

     // warn against using I/O redirection in parallel runs
      if ((inflag == 0) && (universe->nprocs > 1))
        error->warning(FLERR, "Using I/O redirection is unreliable with parallel runs. "
                       "Better to use the -in switch to read input files.");
      utils::flush_buffers(this);
    }

  // universe is one or more worlds, as setup by partition switch
  // split universe communicator into separate world communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    screen = logfile = infile = nullptr;
    if (me == 0) {
      std::string str;
      if (partscreenflag == 0) {
        if (screenflag == 0) {
          str = fmt::format("screen.{}",universe->iworld);
          screen = fopen(str.c_str(),"w");
          if (screen == nullptr)
            error->one(FLERR,"Cannot open screen file {}: {}",str,utils::getsyserror());
        } else if (strcmp(arg[screenflag],"none") == 0) {
          screen = nullptr;
        } else {
          str = fmt::format("{}.{}",arg[screenflag],universe->iworld);
          screen = fopen(str.c_str(),"w");
          if (screen == nullptr)
            error->one(FLERR,"Cannot open screen file {}: {}",arg[screenflag],utils::getsyserror());
        }
      } else if (strcmp(arg[partscreenflag],"none") == 0) {
        screen = nullptr;
      } else {
        str = fmt::format("{}.{}",arg[partscreenflag],universe->iworld);
        screen = fopen(str.c_str(),"w");
        if (screen == nullptr)
          error->one(FLERR,"Cannot open screen file {}: {}",str,utils::getsyserror());
      }

      if (partlogflag == 0) {
        if (logflag == 0) {
          str = fmt::format("log.lammps.{}",universe->iworld);
          logfile = fopen(str.c_str(),"w");
          if (logfile == nullptr)
            error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
        } else if (strcmp(arg[logflag],"none") == 0) {
          logfile = nullptr;
        } else {
          str = fmt::format("{}.{}",arg[logflag],universe->iworld);
          logfile = fopen(str.c_str(),"w");
          if (logfile == nullptr)
            error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
        }
      } else if (strcmp(arg[partlogflag],"none") == 0) {
        logfile = nullptr;
      } else {
        str = fmt::format("{}.{}",arg[partlogflag],universe->iworld);
        logfile = fopen(str.c_str(),"w");
        if (logfile == nullptr)
          error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
      }

      if (strcmp(arg[inflag], "none") != 0) {
        infile = fopen(arg[inflag],"r");
        if (infile == nullptr)
          error->one(FLERR,"Cannot open input script {}: {}",arg[inflag], utils::getsyserror());
      }
    }

    // make all screen and logfile output unbuffered for debugging crashes

    if (nonbufflag) {
      if (universe->uscreen) setbuf(universe->uscreen, nullptr);
      if (universe->ulogfile) setbuf(universe->ulogfile, nullptr);
      if (screen) setbuf(screen, nullptr);
      if (logfile) setbuf(logfile, nullptr);
    }

    // screen and logfile messages for universe and world

    if ((universe->me == 0) && (!helpflag)) {
      const char fmt[] = "LAMMPS ({})\nRunning on {} partitions of processors\n";
      if (universe->uscreen)
        fmt::print(universe->uscreen,fmt,version,universe->nworlds);

      if (universe->ulogfile)
        fmt::print(universe->ulogfile,fmt,version,universe->nworlds);
    }

    if ((me == 0) && (!helpflag))
      utils::logmesg(this,fmt::format("LAMMPS ({})\nProcessor partition = {}\n",
                                      version, universe->iworld));
  }

  // check consistency of datatype settings in lmptype.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in lmptype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR,"Imageint setting in lmptype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in lmptype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in lmptype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_LMP_TAGINT,&mpisize);
  if (mpisize != sizeof(tagint))
      error->all(FLERR,"MPI_LMP_TAGINT and tagint in lmptype.h are not compatible");
  MPI_Type_size(MPI_LMP_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,"MPI_LMP_BIGINT and bigint in lmptype.h are not compatible");

#ifdef LAMMPS_SMALLBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_BIGBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 8 ||
      sizeof(tagint) != 8 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_SMALLSMALL
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 4)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
  // allocate input class now that MPI is fully setup
  input = new Input(this,narg,arg);

  // copy package cmdline arguments

  if (npack > 0) {
    num_package = npack;
    packargs = new char**[npack];
    for (int i=0; i < npack; ++i) {
      int n = plast[i] - pfirst[i];
      packargs[i] = new char*[n+1];
      for (int j=0; j < n; ++j)
        packargs[i][j] = utils::strdup(arg[pfirst[i]+j]);
      packargs[i][n] = nullptr;
    }
    memory->destroy(pfirst);
    memory->destroy(plast);
  }

  // if helpflag set, print help and quit with "success" status
  // otherwise allocate top level classes.

  if (helpflag) {
    error->done(0);
  } else {
    create();
    post_create();
  }

  // if either restart conversion option was used, invoke 2 commands and quit
  // add args between wfirst and wlast to write_data or write_data command
  // add "noinit" to write_data to prevent a system init
  // write_dump will just give a warning message about no init

  if (restart2data || restart2dump) {
    std::string cmd = fmt::format("read_restart {}",restartfile);
    if (restartremap) cmd += " remap\n";
    input->one(cmd);
    if (restart2data) cmd = "write_data ";
    else cmd = "write_dump";
    for (iarg = wfirst; iarg < wlast; iarg++)
       cmd += fmt::format(" {}", arg[iarg]);
    if (restart2data) cmd += " noinit";
    input->one(cmd);
    error->done(0);
  }
}

/** Shut down a LAMMPS simulation instance
 *
 * The LAMMPS destructor shuts down the simulation by deleting top-level class
 * instances, closing screen and log files for the global instance (aka "world")
 * and files and MPI communicators in sub-partitions ("universes"). Then it
 * deletes the fundamental class instances and copies of data inside the class.
 */
LAMMPS::~LAMMPS()
{
  const int me = comm->me;
  destroy();
  if (num_package) {
    for (int i = 0; i < num_package; i++) {
      for (char **ptr = packargs[i]; *ptr != nullptr; ++ptr)
        delete[] *ptr;
      delete[] packargs[i];
    }
    delete[] packargs;
  }
  num_package = 0;
  packargs = nullptr;

  double totalclock = platform::walltime() - initclock;
  if ((me == 0) && (screen || logfile)) {
    int seconds = fmod(totalclock,60.0);
    totalclock  = (totalclock - seconds) / 60.0;
    int minutes = fmod(totalclock,60.0);
    int hours = (totalclock - minutes) / 60.0;
    utils::logmesg(this,fmt::format("Total wall time: {}:{:02d}:{:02d}\n",
                                    hours, minutes, seconds));
  }

  if (universe->nworlds == 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    logfile = nullptr;
    if (screen != stdout) screen = nullptr;
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
    logfile = nullptr;
    if (screen != stdout) screen = nullptr;
  }

  if (infile && infile != stdin) fclose(infile);

  if (world != universe->uworld) MPI_Comm_free(&world);

  // free the MPI comm created by -mpicolor cmdline arg processed in constructor
  // it was passed to universe as if original universe world
  // may have been split later by partitions, universe will free the splits
  // free a copy of uorig here, so check in universe destructor will still work

  MPI_Comm copy = universe->uorig;
  if (external_comm) MPI_Comm_free(&copy);

  delete input;
  delete universe;
  delete error;
  delete memory;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void LAMMPS::create()
{
  force = nullptr;         // Domain->Lattice checks if Force exists

  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()
  comm = new CommBrick(this);
  neighbor = new Neighbor(this);
  domain = new Domain(this);
  atom = new Atom(this);
  atom->create_avec("atomic",0,nullptr,1);

  group = new Group(this);
  force = new Force(this);    // must be after group, to create temperature
  modify = new Modify(this);
  update = new Update(this);  // must be after output, force, neighbor

  // auto-load plugins
#if defined(LMP_PLUGIN)
  plugin_auto_load(this);
#endif
}

/* ----------------------------------------------------------------------
   check suffix consistency with installed packages
   invoke package-specific default package commands
     only invoke if suffix is set and enabled
     also check if suffix2 is set
   called from LAMMPS constructor and after clear() command
     so that package-specific core classes have been instantiated
------------------------------------------------------------------------- */

void LAMMPS::post_create()
{
  if (skiprunflag) input->one("timer timeout 0 every 1");

  // Don't unnecessarily reissue a package command via suffix
  int package_issued = Suffix::NONE;
  // invoke any command-line package commands
  if (num_package) {
    std::string str;
    for (int i = 0; i < num_package; i++) {
      str = "package";
      char *pkg_name = *(packargs[i]);
      if (pkg_name != nullptr) {
        if (strcmp("gpu", pkg_name) == 0) package_issued |= Suffix::GPU;
        if (strcmp("omp", pkg_name) == 0) package_issued |= Suffix::OMP;
        if (strcmp("intel", pkg_name) == 0) package_issued |= Suffix::INTEL;
      }
      for (char **ptr = packargs[i]; *ptr != nullptr; ++ptr) {
        str += " ";
        str += *ptr;
      }
      input->one(str);
    }
  }
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void LAMMPS::init()
{
  update->init();
  force->init();         // pair must come after update due to minimizer
  domain->init();
  atom->init();          // atom must come after force and domain
                         //   atom deletes extra array
                         //   used by fix shear_history::unpack_restart()
                         //     when force->pair->gran_history creates fix
                         //   atom_vec init uses deform_vremap
  modify->init();        // modify must come after update, force, atom, domain
  neighbor->init();      // neighbor must come after force, modify
  comm->init();          // comm must come after force, modify, neighbor, atom
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void LAMMPS::destroy()
{
  // must wipe out all plugins first, if configured
#if defined(LMP_PLUGIN)
  plugin_clear(this);
#endif

  delete update;
  update = nullptr;

  delete neighbor;
  neighbor = nullptr;

  delete force;
  force = nullptr;

  delete group;
  group = nullptr;

  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  modify = nullptr;

  delete comm;            // comm must come after modify
                          //   since fix destructors may access comm
  comm = nullptr;

  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  domain = nullptr;

  delete atom;            // atom must come after modify, neighbor
                          //   since fixes delete callbacks in atom
  atom = nullptr;
  restart_ver = -1;       // reset last restart version id
}

