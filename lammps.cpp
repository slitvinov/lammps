#include "lammps.h"
#include "style_atom.h"
#include "style_command.h"
#include "style_integrate.h"
#include "style_pair.h"
#include "style_region.h"
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
#include "universe.h"
#include "update.h"
#include "variable.h"
#include "version.h"
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
  const char *lmpenv = getenv("LAMMPSHOME");
  if (lmpenv) {
    platform::putenv(fmt::format("PYTHONHOME={}",lmpenv));
  }
#endif
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
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
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
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
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
  if (universe->existflag == 0) universe->add_world(nullptr);
  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions do not match number of allocated processors");
  if (universe->existflag && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");
  if (universe->existflag == 0 && partscreenflag)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");
  if (universe->existflag == 0 && partlogflag)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");
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
      if ((inflag == 0) && (universe->nprocs > 1))
        error->warning(FLERR, "Using I/O redirection is unreliable with parallel runs. "
                       "Better to use the -in switch to read input files.");
      utils::flush_buffers(this);
    }
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
    if (nonbufflag) {
      if (universe->uscreen) setbuf(universe->uscreen, nullptr);
      if (universe->ulogfile) setbuf(universe->ulogfile, nullptr);
      if (screen) setbuf(screen, nullptr);
      if (logfile) setbuf(logfile, nullptr);
    }
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
  input = new Input(this,narg,arg);
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
  if (helpflag) {
    error->done(0);
  } else {
    create();
    post_create();
  }
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
    totalclock = (totalclock - seconds) / 60.0;
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
  MPI_Comm copy = universe->uorig;
  if (external_comm) MPI_Comm_free(&copy);
  delete input;
  delete universe;
  delete error;
  delete memory;
}
void LAMMPS::create()
{
  force = nullptr;
  comm = new CommBrick(this);
  neighbor = new Neighbor(this);
  domain = new Domain(this);
  atom = new Atom(this);
  atom->create_avec("atomic",0,nullptr,1);
  group = new Group(this);
  force = new Force(this);
  modify = new Modify(this);
  update = new Update(this);
}
void LAMMPS::post_create()
{
  if (skiprunflag) input->one("timer timeout 0 every 1");
}
void LAMMPS::init()
{
  update->init();
  force->init();
  domain->init();
  atom->init();
  modify->init();
  neighbor->init();
  comm->init();
}
void LAMMPS::destroy()
{
  delete update;
  update = nullptr;
  delete neighbor;
  neighbor = nullptr;
  delete force;
  force = nullptr;
  delete group;
  group = nullptr;
  delete modify;
  modify = nullptr;
  delete comm;
  comm = nullptr;
  delete domain;
  domain = nullptr;
  delete atom;
  atom = nullptr;
  restart_ver = -1;
}
