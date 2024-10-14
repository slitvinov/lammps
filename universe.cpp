#include "universe.h"
#include "error.h"
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
Universe::~Universe() {
  if (uworld != uorig)
    MPI_Comm_free(&uworld);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(uni2orig);
}
void Universe::reorder(char *style, char *arg) {
  char line[MAXLINE];
  if (uworld != uorig)
    MPI_Comm_free(&uworld);
  if (strcmp(style, "nth") == 0) {
    int n = utils::inumeric(FLERR, arg, false, lmp);
    if (n <= 0)
      error->universe_all(FLERR, "Invalid -reorder N value");
    if (nprocs % n)
      error->universe_all(FLERR, "Nprocs not a multiple of N for -reorder");
    for (int i = 0; i < nprocs; i++) {
      if (i < (n - 1) * nprocs / n)
        uni2orig[i] = i / (n - 1) * n + (i % (n - 1));
      else
        uni2orig[i] = (i - (n - 1) * nprocs / n) * n + n - 1;
    }
  } else if (strcmp(style, "custom") == 0) {
    if (me == 0) {
      FILE *fp = fopen(arg, "r");
      if (fp == nullptr)
        error->universe_one(FLERR, fmt::format("Cannot open -reorder "
                                               "file {}: {}",
                                               arg, utils::getsyserror()));
      char *ptr;
      if (!fgets(line, MAXLINE, fp))
        error->one(FLERR, "Unexpected end of -reorder file");
      while (true) {
        if ((ptr = strchr(line, '#')))
          *ptr = '\0';
        if (strspn(line, " \t\n\r") != strlen(line))
          break;
        if (!fgets(line, MAXLINE, fp))
          error->one(FLERR, "Unexpected end of -reorder file");
      }
      int me_orig, me_new, rv;
      rv = sscanf(line, "%d %d", &me_orig, &me_new);
      if (me_orig < 0 || me_orig >= nprocs || me_new < 0 || me_new >= nprocs ||
          rv != 2)
        error->one(FLERR,
                   "Invalid entry '{} {}' in -reorder "
                   "file",
                   me_orig, me_new);
      uni2orig[me_new] = me_orig;
      for (int i = 1; i < nprocs; i++) {
        if (!fgets(line, MAXLINE, fp))
          error->one(FLERR, "Unexpected end of -reorder file");
        rv = sscanf(line, "%d %d", &me_orig, &me_new);
        if (me_orig < 0 || me_orig >= nprocs || me_new < 0 ||
            me_new >= nprocs || rv != 2)
          error->one(FLERR,
                     "Invalid entry '{} {}' in -reorder "
                     "file",
                     me_orig, me_new);
        uni2orig[me_new] = me_orig;
      }
      fclose(fp);
    }
    MPI_Bcast(uni2orig, nprocs, MPI_INT, 0, uorig);
  } else
    error->universe_all(FLERR, "Invalid command-line argument");
  int ome, key;
  MPI_Comm_rank(uorig, &ome);
  for (int i = 0; i < nprocs; i++)
    if (uni2orig[i] == ome)
      key = i;
  MPI_Comm_split(uorig, 0, key, &uworld);
  MPI_Comm_rank(uworld, &me);
  MPI_Comm_size(uworld, &nprocs);
}
void Universe::add_world(char *str) {
  int n, nper;
  n = 1;
  nper = 0;
  if (str != nullptr) {
    bool valid = true;
    std::string part(str);
    if (part.size() == 0)
      valid = false;
    if (part.find_first_not_of("0123456789x") != std::string::npos)
      valid = false;
    if (valid) {
      std::size_t found = part.find_first_of('x');
      if ((found == 0) || (found == (part.size() - 1))) {
        valid = false;
      } else if (found == std::string::npos) {
        nper = atoi(part.c_str());
      } else {
        n = atoi(part.substr(0, found).c_str());
        nper = atoi(part.substr(found + 1).c_str());
      }
    }
    if (n < 1 || nper < 1)
      valid = false;
    if (!valid)
      error->universe_all(FLERR,
                          fmt::format("Invalid partition string '{}'", str));
  } else
    nper = nprocs;
  memory->grow(procs_per_world, nworlds + n, "universe:procs_per_world");
  memory->grow(root_proc, (nworlds + n), "universe:root_proc");
  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0)
      root_proc[nworlds] = 0;
    else
      root_proc[nworlds] =
          root_proc[nworlds - 1] + procs_per_world[nworlds - 1];
    if (me >= root_proc[nworlds])
      iworld = nworlds;
    nworlds++;
  }
}
int Universe::consistent() {
  int n = 0;
  for (int i = 0; i < nworlds; i++)
    n += procs_per_world[i];
  if (n == nprocs)
    return 1;
  else
    return 0;
}
