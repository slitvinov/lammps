#include "input.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "comm_brick.h"
#include "command.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "integrate.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "style_command.h"
#include "universe.h"
#include "update.h"
#include <cctype>
#include <cerrno>
#include <cstring>
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
  echo_screen = 0;
  echo_log = 1;
  label_active = 0;
  labelstr = nullptr;
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
#define COMMAND_CLASS
#define CommandStyle(key, Class) (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
}
Input::~Input() {
  memory->sfree(line);
  memory->sfree(copy);
  memory->sfree(work);
  delete[] labelstr;
  memory->sfree(arg);
  delete[] infiles;
  delete command_map;
}
void Input::file() {
  int m, n, mstart, ntriple, endfile;
  while (true) {
    if (me == 0) {
      ntriple = 0;
      endfile = 0;
      m = 0;
      while (true) {
        if (infile == nullptr) {
          n = 0;
          break;
        }
        mstart = m;
        while (true) {
          if (maxline - m < 2)
            reallocate(line, maxline, 0);
          if (fgets(&line[m], maxline - m, infile) == nullptr) {
            endfile = 1;
            if (m)
              n = strlen(line) + 1;
            else
              n = 0;
            break;
          }
          m += strlen(&line[m]);
          if (line[m - 1] != '\n')
            continue;
          break;
        }
        if (endfile)
          break;
        ntriple += numtriple(&line[mstart]);
        m--;
        while (m >= 0 && isspace(line[m]))
          m--;
        if (m >= 0 && line[m] == '&')
          continue;
        if (ntriple % 2) {
          line[m + 1] = '\n';
          m += 2;
          continue;
        }
        line[m + 1] = '\0';
        n = m + 2;
        break;
      }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n == 0) {
      if (label_active)
        error->all(FLERR, "Label wasn't found in input script");
      break;
    }
    if (n > maxline)
      reallocate(line, maxline, n);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);
    if (me == 0 && label_active == 0) {
      if (echo_screen && screen)
        fprintf(screen, "%s\n", line);
      if (echo_log && logfile)
        fprintf(logfile, "%s\n", line);
    }
    parse();
    if (command == nullptr)
      continue;
    if (label_active && strcmp(command, "label") != 0)
      continue;
    if (execute_command() && line)
      error->all(FLERR, "Unknown command: {}", line);
  }
}
void Input::write_echo(const std::string &txt) {
  if (me == 0) {
    if (echo_screen && screen)
      fputs(txt.c_str(), screen);
    if (echo_log && logfile)
      fputs(txt.c_str(), logfile);
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
  if (!label_active)
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
  if (mycmd == "log")
    log();
  else if (mycmd == "comm_modify")
    comm_modify();
  else if (mycmd == "fix")
    fix();
  else if (mycmd == "fix_modify")
    fix_modify();
  else if (mycmd == "lattice")
    lattice();
  else if (mycmd == "mass")
    mass();
  else if (mycmd == "neigh_modify")
    neigh_modify();
  else if (mycmd == "neighbor")
    neighbor_command();
  else if (mycmd == "pair_coeff")
    pair_coeff();
  else if (mycmd == "pair_style")
    pair_style();
  else if (mycmd == "region")
    region();
  else if (mycmd == "timestep")
    timestep();
  else
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
void Input::log() {
  if ((narg < 1) || (narg > 2))
    error->all(FLERR,
               "Illegal log command: expected 1 or 2 argument(s) but found {}",
               narg);
  int appendflag = 0;
  if (narg == 2) {
    if (strcmp(arg[1], "append") == 0)
      appendflag = 1;
    else
      error->all(FLERR, "Unknown log keyword: {}", arg[1]);
  }
  if (me == 0) {
    if (logfile)
      fclose(logfile);
    if (strcmp(arg[0], "none") == 0)
      logfile = nullptr;
    else {
      if (appendflag)
        logfile = fopen(arg[0], "a");
      else
        logfile = fopen(arg[0], "w");
      if (logfile == nullptr)
        error->one(FLERR, "Cannot open logfile {}: {}", arg[0],
                   utils::getsyserror());
    }
    if (universe->nworlds == 1)
      universe->ulogfile = logfile;
  }
}
void Input::comm_modify() { comm->modify_params(narg, arg); }
void Input::fix() { modify->add_fix(narg, arg); }
void Input::fix_modify() { modify->modify_fix(narg, arg); }
void Input::lattice() { domain->set_lattice(narg, arg); }
void Input::mass() {
  if (narg != 2)
    error->all(FLERR, "Illegal mass command: expected 2 arguments but found {}",
               narg);
  if (domain->box_exist == 0)
    error->all(FLERR, "Mass command before simulation box is defined");
  atom->set_mass(FLERR, narg, arg);
}
void Input::neigh_modify() { neighbor->modify_params(narg, arg); }
void Input::neighbor_command() { neighbor->set(narg, arg); }
void Input::pair_coeff() {
  if (domain->box_exist == 0)
    error->all(FLERR, "Pair_coeff command before simulation box is defined");
  if (force->pair == nullptr)
    error->all(FLERR, "Pair_coeff command without a pair style");
  if (narg < 2)
    utils::missing_cmd_args(FLERR, "pair_coeff", error);
  if (force->pair->one_coeff &&
      ((strcmp(arg[0], "*") != 0) || (strcmp(arg[1], "*") != 0)))
    error->all(FLERR, "Pair_coeff must start with * * for pair style {}",
               force->pair_style);
  int itype, jtype;
  if (utils::strmatch(arg[0], "^\\d+$") && utils::strmatch(arg[1], "^\\d+$")) {
    itype = utils::inumeric(FLERR, arg[0], false, lmp);
    jtype = utils::inumeric(FLERR, arg[1], false, lmp);
    if (jtype < itype) {
      char *str = arg[0];
      arg[0] = arg[1];
      arg[1] = str;
    }
  }
  force->pair->coeff(narg, arg);
}
void Input::pair_style() {
  if (narg < 1)
    utils::missing_cmd_args(FLERR, "pair_style", error);
  if (force->pair) {
    std::string style = arg[0];
    int match = 0;
    if (style == force->pair_style)
      match = 1;
    if (match) {
      force->pair->settings(narg - 1, &arg[1]);
      return;
    }
  }
  force->create_pair(arg[0], 1);
  if (force->pair)
    force->pair->settings(narg - 1, &arg[1]);
}
void Input::region() { domain->add_region(narg, arg); }
void Input::timestep() {
  if (narg != 1)
    error->all(FLERR, "Illegal timestep command");
  update->update_time();
  update->dt = utils::numeric(FLERR, arg[0], false, lmp);
  update->dt_default = 0;
  if (update->first_update == 0)
    return;
  int respaflag = 0;
  if (utils::strmatch(update->integrate_style, "^respa"))
    respaflag = 1;
  if (respaflag)
    update->integrate->reset_dt();
  if (force->pair)
    force->pair->reset_dt();
  for (auto &ifix : modify->get_fix_list())
    ifix->reset_dt();
}
