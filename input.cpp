#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <string>
#include <cmath>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "input.h"
#include "atom.h"
#include "atom_vec.h"
#include "pointers.h"
#include "comm.h"
#include "comm_brick.h"
#include "pointers.h"
#include "command.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "create_atoms.h"
#include "create_box.h"
#include "run.h"
#include "universe.h"
#include "update.h"
using namespace LAMMPS_NS;
#define DELTALINE 256
#define DELTA 4
template <typename T> static Command *command_creator(LAMMPS *lmp) {
  return new T(lmp);
}
Input::Input(LAMMPS *lmp, int argc, char **argv) : Pointers(lmp) {
  MPI_Comm_rank(world, &me);
  maxline = maxcopy = maxwork = 0;
  line = copy = work = nullptr;
  narg = maxarg = 0;
  arg = nullptr;
  jump_skip = 0;
  utf8_warn = true;
  if (me == 0) {
    nfile = 1;
    maxfile = 16;
    infiles = new FILE *[maxfile];
    infiles[0] = infile;
  } else
    infiles = nullptr;
  command_map = new CommandCreatorMap();
  (*command_map)["run"] = &command_creator<Run>;
  (*command_map)["create_atoms"] = &command_creator<CreateAtoms>;
  (*command_map)["create_box"] = &command_creator<CreateBox>;
}
void Input::file() {
  int m, n, mstart, ntriple, endfile;
  while (true) {
    if (me == 0) {
      ntriple = 0;
      endfile = 0;
      m = 0;
      while (true) {
	mstart = m;
	while (true) {
	  if (maxline - m < 2)
	    reallocate(line, maxline, 0);
	  if (fgets(&line[m], maxline - m, infile) == nullptr) {
	    endfile = 1;
	    n = 0;
	    break;
	  }
	  m += strlen(&line[m]);
	  break;
	}
	if (endfile)
	  break;
	ntriple += numtriple(&line[mstart]);
	m--;
	while (m >= 0 && isspace(line[m]))
	  m--;
	line[m + 1] = '\0';
	n = m + 2;
	break;
      }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n == 0) {
      break;
    }
    if (n > maxline)
      reallocate(line, maxline, n);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);
    parse();
    if (command == nullptr)
      continue;
    execute_command();
  }
}
void Input::parse() {
  int n = strlen(line) + 1;
  if (n > maxcopy)
    reallocate(copy, maxcopy, n);
  strcpy(copy, line);
  char *ptrmatch;
  char *ptr = copy;
  while (*ptr) {
    ptr++;
  }
  substitute(copy, work, maxcopy, maxwork, 1);
  char *next;
  command = nextword(copy, &next);
  if (command == nullptr)
    return;
  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg =
	  (char **)memory->srealloc(arg, maxarg * sizeof(char *), "input:arg");
    }
    arg[narg] = nextword(ptr, &next);
    if (!arg[narg])
      break;
    narg++;
    ptr = next;
  }
}
char *Input::nextword(char *str, char **next) {
  char *start, *stop;
  start = &str[strspn(str, " \t\n\v\f\r")];
  if (*start == '\0')
    return nullptr;
  stop = &start[strcspn(start, " \t\n\v\f\r")];
  if (*stop == '\0')
    *next = stop;
  else
    *next = stop + 1;
  *stop = '\0';
  return start;
}
void Input::substitute(char *&str, char *&str2, int &max, int &max2, int flag) {
  int i, n, paren_count, nchars;
  ;
  char immediate[256];
  char *var, *value, *beyond;
  int quoteflag = 0;
  char *ptrmatch;
  char *ptr = str;
  n = strlen(str) + 1;
  if (n > max2)
    reallocate(str2, max2, n);
  *str2 = '\0';
  char *ptr2 = str2;
  while (*ptr) {
    *ptr2++ = *ptr++;
    *ptr2 = '\0';
  }
  if (max2 > max)
    reallocate(str, max, max2);
  strcpy(str, str2);
}
int Input::numtriple(char *line) {
  int count = 0;
  char *ptr = line;
  return count;
}
void Input::reallocate(char *&str, int &max, int n) {
  if (n) {
    while (n > max)
      max += DELTALINE;
  } else
    max += DELTALINE;
  str = (char *)memory->srealloc(str, max * sizeof(char), "input:str");
}
int Input::execute_command() {
  int flag = 1;
  std::string mycmd = command;
  if (mycmd == "comm_modify")
    comm->modify_params(narg, arg);
  else if (mycmd == "fix")
    modify->add_fix(narg, arg);
  else if (mycmd == "mass")
    atom->set_mass(FLERR, narg, arg);
  else if (mycmd == "neigh_modify")
    neighbor->modify_params(narg, arg);
  else if (mycmd == "neighbor")
    neighbor->set(narg, arg);
  else if (mycmd == "pair_coeff") {
    force->pair->coeff(narg, arg);
  } else if (mycmd == "pair_style") {
    force->create_pair(arg[0], 1);
    if (force->pair)
      force->pair->settings(narg - 1, &arg[1]);
  } else if (mycmd == "region")
    domain->add_region(narg, arg);
  else if (mycmd == "timestep") {
    update->update_time();
    update->dt = utils::numeric(FLERR, arg[0], false, lmp);
    update->dt_default = 0;
  } else
    flag = 0;
  if (flag)
    return 0;
  if (command_map->find(mycmd) != command_map->end()) {
    CommandCreator &command_creator = (*command_map)[mycmd];
    Command *cmd = command_creator(lmp);
    cmd->command(narg, arg);
    delete cmd;
    return 0;
  }
  return -1;
}
