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
#include "variable.h"
#include <cstring>
#include <cerrno>
#include <cctype>
using namespace LAMMPS_NS;
#define DELTALINE 256
#define DELTA 4
template <typename T> static Command *command_creator(LAMMPS *lmp)
{
  return new T(lmp);
}
Input::Input(LAMMPS *lmp, int argc, char **argv) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
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
  } else infiles = nullptr;
  variable = new Variable(lmp);
  command_map = new CommandCreatorMap();
#define COMMAND_CLASS 
#define CommandStyle(key,Class) \
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
}
Input::~Input()
{
  memory->sfree(line);
  memory->sfree(copy);
  memory->sfree(work);
  delete[] labelstr;
  memory->sfree(arg);
  delete[] infiles;
  delete variable;
  delete command_map;
}
void Input::file()
{
  int m,n,mstart,ntriple,endfile;
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
          if (maxline-m < 2) reallocate(line,maxline,0);
          if (fgets(&line[m],maxline-m,infile) == nullptr) {
            endfile = 1;
            if (m) n = strlen(line) + 1;
            else n = 0;
            break;
          }
          m += strlen(&line[m]);
          if (line[m-1] != '\n') continue;
          break;
        }
        if (endfile) break;
        ntriple += numtriple(&line[mstart]);
        m--;
        while (m >= 0 && isspace(line[m])) m--;
        if (m >= 0 && line[m] == '&') continue;
        if (ntriple % 2) {
          line[m+1] = '\n';
          m += 2;
          continue;
        }
        line[m+1] = '\0';
        n = m+2;
        break;
      }
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all(FLERR,"Label wasn't found in input script");
      break;
    }
    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s\n",line);
      if (echo_log && logfile) fprintf(logfile,"%s\n",line);
    }
    parse();
    if (command == nullptr) continue;
    if (label_active && strcmp(command,"label") != 0) continue;
    if (execute_command() && line)
      error->all(FLERR,"Unknown command: {}",line);
  }
}
void Input::write_echo(const std::string &txt)
{
  if (me == 0) {
    if (echo_screen && screen) fputs(txt.c_str(),screen);
    if (echo_log && logfile) fputs(txt.c_str(),logfile);
  }
}
void Input::parse()
{
  int n = strlen(line) + 1;
  if (n > maxcopy) reallocate(copy,maxcopy,n);
  strcpy(copy,line);
  char *ptrmatch;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#') {
      *ptr = '\0';
      break;
    }
    if (*ptr == '\'') {
      ptrmatch = strchr(ptr+1,'\'');
      if (ptrmatch == nullptr)
        error->all(FLERR,"Unmatched single quote in command");
      ptr = ptrmatch + 1;
    } else if (*ptr == '"') {
      if (strstr(ptr,"\"\"\"") == ptr) {
        ptrmatch = strstr(ptr+3,"\"\"\"");
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched triple quote in command");
        ptr = ptrmatch + 3;
      } else {
        ptrmatch = strchr(ptr+1,'"');
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched double quote in command");
        ptr = ptrmatch + 1;
      }
    } else ptr++;
  }
  if (utils::has_utf8(copy)) {
    std::string buf = utils::utf8_subst(copy);
    strcpy(copy,buf.c_str());
    if (utf8_warn && (comm->me == 0))
      error->warning(FLERR,"Detected non-ASCII characters in input. "
                     "Will try to continue by replacing with ASCII "
                     "equivalents where known.");
    utf8_warn = false;
  }
  if (!label_active) substitute(copy,work,maxcopy,maxwork,1);
  char *next;
  command = nextword(copy,&next);
  if (command == nullptr) return;
  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}
char *Input::nextword(char *str, char **next)
{
  char *start,*stop;
  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return nullptr;
  if (strstr(start,"\"\"\"") == start) {
    stop = strstr(&start[3],"\"\"\"");
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start += 3;
    *next = stop+3;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by white-space");
  } else if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start++;
    *next = stop+1;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by white-space");
  } else {
    stop = &start[strcspn(start," \t\n\v\f\r")];
    if (*stop == '\0') *next = stop;
    else *next = stop+1;
  }
  *stop = '\0';
  return start;
}
void Input::substitute(char *&str, char *&str2, int &max, int &max2, int flag)
{
  int i,n,paren_count,nchars;;
  char immediate[256];
  char *var,*value,*beyond;
  int quoteflag = 0;
  char *ptrmatch;
  char *ptr = str;
  n = strlen(str) + 1;
  if (n > max2) reallocate(str2,max2,n);
  *str2 = '\0';
  char *ptr2 = str2;
  while (*ptr) {
    if (*ptr == '$' && !quoteflag) {
      if (*(ptr+1) == '{') {
        var = ptr+2;
        i = 0;
        while (var[i] != '\0' && var[i] != '}') i++;
        if (var[i] == '\0') error->one(FLERR,"Invalid variable name");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
        value = variable->retrieve(var);
      } else if (*(ptr+1) == '(') {
        var = ptr+2;
        paren_count = 0;
        i = 0;
        while (var[i] != '\0' && (var[i] != ')' || paren_count != 0)) {
          switch (var[i]) {
          case '(': paren_count++; break;
          case ')': paren_count--; break;
          default: ;
          }
          i++;
        }
        if (var[i] == '\0') error->one(FLERR,"Invalid immediate variable");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
        char fmtstr[64] = "%.20g";
        char *fmtflag;
        if ((fmtflag=strrchr(var, ':')) && (fmtflag[1]=='%')) {
          strncpy(fmtstr,&fmtflag[1],sizeof(fmtstr)-1);
          *fmtflag='\0';
        }
        if (!utils::strmatch(fmtstr,"%[0-9 ]*\\.[0-9]+[efgEFG]"))
          error->all(FLERR,"Incorrect conversion in format string");
        snprintf(immediate,256,fmtstr,variable->compute_equal(var));
        value = immediate;
      } else {
        var = ptr;
        var[0] = var[1];
        var[1] = '\0';
        beyond = ptr + 2;
        value = variable->retrieve(var);
      }
      if (value == nullptr)
        error->one(FLERR,"Substitution for illegal variable {}",var);
      n = strlen(str2) + strlen(value) + strlen(beyond) + 1;
      if (n > max2) reallocate(str2,max2,n);
      strcat(str2,value);
      ptr2 = str2 + strlen(str2);
      ptr = beyond;
      if (flag && me == 0 && label_active == 0) {
        if (echo_screen && screen) fprintf(screen,"%s%s\n",str2,beyond);
        if (echo_log && logfile) fprintf(logfile,"%s%s\n",str2,beyond);
      }
    } else if (*ptr == '\'') {
      ptrmatch = strchr(ptr+1,'\'');
      if (ptrmatch == nullptr)
        error->all(FLERR,"Unmatched single quote in command");
      nchars = ptrmatch+1 - ptr;
      strncpy(ptr2,ptr,nchars);
      ptr += nchars;
      ptr2 += nchars;
    } else if (*ptr == '"') {
      if (strstr(ptr,"\"\"\"") == ptr) {
        ptrmatch = strstr(ptr+3,"\"\"\"");
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched triple quote in command");
        nchars = ptrmatch+3 - ptr;
        strncpy(ptr2,ptr,nchars);
        ptr += nchars;
        ptr2 += nchars;
      } else {
        ptrmatch = strchr(ptr+1,'"');
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched double quote in command");
        nchars = ptrmatch+1 - ptr;
        strncpy(ptr2,ptr,nchars);
        ptr += nchars;
        ptr2 += nchars;
      }
    } else *ptr2++ = *ptr++;
    *ptr2 = '\0';
  }
  if (max2 > max) reallocate(str,max,max2);
  strcpy(str,str2);
}
int Input::numtriple(char *line)
{
  int count = 0;
  char *ptr = line;
  while ((ptr = strstr(ptr,"\"\"\""))) {
    ptr += 3;
    count++;
  }
  return count;
}
void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;
  str = (char *) memory->srealloc(str,max*sizeof(char),"input:str");
}
int Input::execute_command()
{
  int flag = 1;
  std::string mycmd = command;
  if (mycmd == "log") log();
  else if (mycmd == "print") print();
  else if (mycmd == "atom_modify") atom_modify();
  else if (mycmd == "atom_style") atom_style();
  else if (mycmd == "boundary") boundary();
  else if (mycmd == "comm_modify") comm_modify();
  else if (mycmd == "comm_style") comm_style();
  else if (mycmd == "compute") compute();
  else if (mycmd == "compute_modify") compute_modify();
  else if (mycmd == "dielectric") dielectric();
  else if (mycmd == "dimension") dimension();
  else if (mycmd == "fix") fix();
  else if (mycmd == "fix_modify") fix_modify();
  else if (mycmd == "group") group_command();
  else if (mycmd == "lattice") lattice();
  else if (mycmd == "mass") mass();
  else if (mycmd == "neigh_modify") neigh_modify();
  else if (mycmd == "neighbor") neighbor_command();
  else if (mycmd == "newton") newton();
  else if (mycmd == "pair_coeff") pair_coeff();
  else if (mycmd == "pair_modify") pair_modify();
  else if (mycmd == "pair_style") pair_style();
  else if (mycmd == "pair_write") pair_write();
  else if (mycmd == "processors") processors();
  else if (mycmd == "region") region();
  else if (mycmd == "run_style") run_style();
  else if (mycmd == "timestep") timestep();
  else flag = 0;
  if (flag) return 0;
  if (mycmd == "reset_atoms") flag = meta(mycmd);
  if (flag) return 0;
  if (command_map->find(mycmd) != command_map->end()) {
    CommandCreator &command_creator = (*command_map)[mycmd];
    Command *cmd = command_creator(lmp);
    cmd->command(narg,arg);
    delete cmd;
    return 0;
  }
  return -1;
}
void Input::log()
{
  if ((narg < 1) || (narg > 2)) error->all(FLERR,"Illegal log command: expected 1 or 2 argument(s) but found {}", narg);
  int appendflag = 0;
  if (narg == 2) {
    if (strcmp(arg[1],"append") == 0) appendflag = 1;
    else error->all(FLERR,"Unknown log keyword: {}", arg[1]);
  }
  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = nullptr;
    else {
      if (appendflag) logfile = fopen(arg[0],"a");
      else logfile = fopen(arg[0],"w");
      if (logfile == nullptr)
        error->one(FLERR,"Cannot open logfile {}: {}",
                                     arg[0], utils::getsyserror());
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}
void Input::print()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "print", error);
  int n = strlen(arg[0]) + 1;
  if (n > maxline) reallocate(line,maxline,n);
  strcpy(line,arg[0]);
  substitute(line,work,maxline,maxwork,0);
  FILE *fp = nullptr;
  int screenflag = 1;
  int universeflag = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal print {} command: missing argument(s)", arg[iarg]);
      if (me == 0) {
        if (fp != nullptr) fclose(fp);
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open print file {}: {}", arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "print screen", error);
      screenflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"universe") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "print universe", error);
      universeflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown print keyword: {}", arg[iarg]);
  }
  if (me == 0) {
    if (screenflag && screen) fprintf(screen,"%s\n",line);
    if (screenflag && logfile) fprintf(logfile,"%s\n",line);
    if (fp) {
      fprintf(fp,"%s\n",line);
      fclose(fp);
    }
  }
  if (universeflag && (universe->me == 0)) {
    if (universe->uscreen) fprintf(universe->uscreen, "%s\n",line);
    if (universe->ulogfile) fprintf(universe->ulogfile,"%s\n",line);
  }
}
void Input::atom_modify()
{
  atom->modify_params(narg,arg);
}
void Input::atom_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "atom_style", error);
  if (domain->box_exist)
    error->all(FLERR,"Atom_style command after simulation box is defined");
  atom->create_avec(arg[0],narg-1,&arg[1],1);
}
void Input::boundary()
{
  if (domain->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");
  domain->set_boundary(narg,arg,0);
}
void Input::comm_modify()
{
  comm->modify_params(narg,arg);
}
void Input::comm_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "comm_style", error);
  if (strcmp(arg[0],"brick") == 0) {
    if (comm->style == Comm::BRICK) return;
    Comm *oldcomm = comm;
    comm = new CommBrick(lmp,oldcomm);
    delete oldcomm;
  } else error->all(FLERR,"Unknown comm_style argument: {}", arg[0]);
}
void Input::compute()
{
  modify->add_compute(narg,arg);
}
void Input::compute_modify()
{
  modify->modify_compute(narg,arg);
}
void Input::dielectric()
{
  if (narg != 1) error->all(FLERR,"Illegal dielectric command");
  force->dielectric = utils::numeric(FLERR,arg[0],false,lmp);
}
void Input::dimension()
{
  if (narg != 1) error->all(FLERR, "Dimension command expects exactly 1 argument");
  if (domain->box_exist)
    error->all(FLERR,"Dimension command after simulation box is defined");
  domain->dimension = utils::inumeric(FLERR,arg[0],false,lmp);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR, "Invalid dimension argument: {}", arg[0]);
  for (auto &c : modify->get_compute_list()) c->reset_extra_dof();
}
void Input::fix()
{
  modify->add_fix(narg,arg);
}
void Input::fix_modify()
{
  modify->modify_fix(narg,arg);
}
void Input::group_command()
{
  group->assign(narg,arg);
}
void Input::lattice()
{
  domain->set_lattice(narg,arg);
}
void Input::mass()
{
  if (narg != 2) error->all(FLERR,"Illegal mass command: expected 2 arguments but found {}", narg);
  if (domain->box_exist == 0)
    error->all(FLERR,"Mass command before simulation box is defined");
  atom->set_mass(FLERR,narg,arg);
}
void Input::neigh_modify()
{
  neighbor->modify_params(narg,arg);
}
void Input::neighbor_command()
{
  neighbor->set(narg,arg);
}
void Input::newton()
{
  int newton_pair=1,newton_bond=1;
  if (narg == 1) {
    newton_pair = newton_bond = utils::logical(FLERR,arg[0],false,lmp);
  } else if (narg == 2) {
    newton_pair = utils::logical(FLERR,arg[0],false,lmp);
    newton_bond = utils::logical(FLERR,arg[1],false,lmp);
  } else error->all(FLERR,"Illegal newton command");
  force->newton_pair = newton_pair;
  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;
}
void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Pair_coeff command before simulation box is defined");
  if (force->pair == nullptr) error->all(FLERR,"Pair_coeff command without a pair style");
  if (narg < 2) utils::missing_cmd_args(FLERR,"pair_coeff", error);
  if (force->pair->one_coeff && ((strcmp(arg[0],"*") != 0) || (strcmp(arg[1],"*") != 0)))
    error->all(FLERR,"Pair_coeff must start with * * for pair style {}", force->pair_style);
  int itype,jtype;
  if (utils::strmatch(arg[0],"^\\d+$") && utils::strmatch(arg[1],"^\\d+$")) {
    itype = utils::inumeric(FLERR,arg[0],false,lmp);
    jtype = utils::inumeric(FLERR,arg[1],false,lmp);
    if (jtype < itype) {
      char *str = arg[0];
      arg[0] = arg[1];
      arg[1] = str;
    }
  }
  force->pair->coeff(narg,arg);
}
void Input::pair_modify()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Pair_modify command before pair_style is defined");
  force->pair->modify_params(narg,arg);
}
void Input::pair_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "pair_style", error);
  if (force->pair) {
    std::string style = arg[0];
    int match = 0;
    if (style == force->pair_style) match = 1;
    if (match) {
      force->pair->settings(narg-1,&arg[1]);
      return;
    }
  }
  force->create_pair(arg[0],1);
  if (force->pair) force->pair->settings(narg-1,&arg[1]);
}
void Input::pair_write()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Pair_write command before pair_style is defined");
  force->pair->write_file(narg,arg);
}
void Input::processors()
{
  if (domain->box_exist)
    error->all(FLERR,"Processors command after simulation box is defined");
  comm->set_processors(narg,arg);
}
void Input::region()
{
  domain->add_region(narg,arg);
}
void Input::run_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Run_style command before simulation box is defined");
  update->create_integrate(narg,arg,1);
}
void Input::timestep()
{
  if (narg != 1) error->all(FLERR,"Illegal timestep command");
  update->update_time();
  update->dt = utils::numeric(FLERR,arg[0],false,lmp);
  update->dt_default = 0;
  if (update->first_update == 0) return;
  int respaflag = 0;
  if (utils::strmatch(update->integrate_style, "^respa")) respaflag = 1;
  if (respaflag) update->integrate->reset_dt();
  if (force->pair) force->pair->reset_dt();
  for (auto &ifix : modify->get_fix_list()) ifix->reset_dt();
}
int Input::meta(const std::string &prefix)
{
  auto mycmd = fmt::format("{}_{}", utils::uppercase(prefix), utils::uppercase(arg[0]));
  if (command_map->find(mycmd) != command_map->end()) {
    CommandCreator &command_creator = (*command_map)[mycmd];
    Command *cmd = command_creator(lmp);
    cmd->command(narg-1,arg+1);
    delete cmd;
    return 1;
  } else return 0;
}
