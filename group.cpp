#include "group.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "tokenizer.h"
#include "variable.h"
#include "exceptions.h"
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
using namespace LAMMPS_NS;
static constexpr int MAX_GROUP = 32;
static constexpr double EPSILON = 1.0e-6;
enum{NONE,TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};
#define BIG 1.0e20
Group::Group(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  names = new char*[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];
  for (int i = 0; i < MAX_GROUP; i++) names[i] = nullptr;
  for (int i = 0; i < MAX_GROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++) inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++) dynamic[i] = 0;
  names[0] = utils::strdup("all");
  ngroup = 1;
}
Group::~Group()
{
  for (int i = 0; i < MAX_GROUP; i++) delete [] names[i];
  delete [] names;
  delete [] bitmask;
  delete [] inversemask;
  delete [] dynamic;
}
void Group::assign(int narg, char **arg)
{
  int i;
  if (domain->box_exist == 0)
    error->all(FLERR,"Group command before simulation box is defined");
  if (narg < 2) utils::missing_cmd_args(FLERR, "group", error);
  if (strcmp(arg[1],"delete") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal group delete command: too many arguments");
    int igroup = find(arg[0]);
    if (igroup == -1) error->all(FLERR,"Could not find group delete group ID {}",arg[0]);
    if (igroup == 0) error->all(FLERR,"Cannot delete group all");
    for (const auto &i : modify->get_fix_list())
      if (i->igroup == igroup)
        error->all(FLERR,"Cannot delete group {} currently used by fix ID {}",arg[0],i->id);
    for (const auto &i : modify->get_compute_list())
      if (i->igroup == igroup)
        error->all(FLERR,"Cannot delete group {} currently used by compute ID {}",arg[0],i->id);
    if (atom->firstgroupname && strcmp(arg[0],atom->firstgroupname) == 0)
      error->all(FLERR,"Cannot delete group {} currently used by atom_modify first",arg[0]);
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int bits = inversemask[igroup];
    for (i = 0; i < nlocal; i++) mask[i] &= bits;
    if (dynamic[igroup])
      modify->delete_fix(std::string("GROUP_") + names[igroup]);
    delete [] names[igroup];
    names[igroup] = nullptr;
    dynamic[igroup] = 0;
    ngroup--;
    return;
  }
  if (strcmp(arg[1],"clear") == 0) {
    int igroup = find(arg[0]);
    if (igroup == -1) error->all (FLERR,"Could not find group clear group ID {}",arg[0]);
    if (igroup == 0) error->all (FLERR,"Cannot clear group all");
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int bits = inversemask[igroup];
    for (i = 0; i < nlocal; i++) mask[i] &= bits;
    return;
  }
  int igroup = find(arg[0]);
  bool created = false;
  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR,"Too many groups (max {})", MAX_GROUP);
    igroup = find_unused();
    names[igroup] = utils::strdup(arg[0]);
    ngroup++;
    created = true;
  }
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int bit = bitmask[igroup];
  try {
    if (strcmp(arg[1],"region") == 0) {
      if (narg != 3) error->all(FLERR,"Illegal group region command");
      auto region = domain->get_region_by_id(arg[2]);
      if (!region) error->all(FLERR,"Region {} for group region does not exist", arg[2]);
      region->init();
      region->prematch();
      for (i = 0; i < nlocal; i++)
        if (region->match(x[i][0],x[i][1],x[i][2]))
          mask[i] |= bit;
    } else if (strcmp(arg[1],"empty") == 0) {
      if (narg != 2) error->all(FLERR,"Illegal group empty command");
    } else if (strcmp(arg[1],"type") == 0 || strcmp(arg[1],"molecule") == 0 ||
               strcmp(arg[1],"id") == 0) {
      if (narg < 3) utils::missing_cmd_args(FLERR, std::string("group ")+arg[1], error);
      int category=NONE;
      if (strcmp(arg[1],"type") == 0) category = TYPE;
      else if (strcmp(arg[1],"molecule") == 0) category = MOLECULE;
      else if (strcmp(arg[1],"id") == 0) category = ID;
      if ((category == ID) && (!atom->tag_enable))
        error->all(FLERR,"Group id command requires atom IDs");
      if (narg > 3 &&
          (strcmp(arg[2],"<") == 0 || strcmp(arg[2],">") == 0 ||
           strcmp(arg[2],"<=") == 0 || strcmp(arg[2],">=") == 0 ||
           strcmp(arg[2],"==") == 0 || strcmp(arg[2],"!=") == 0 ||
           strcmp(arg[2],"<>") == 0)) {
        int condition = -1;
        if (strcmp(arg[2],"<") == 0) condition = LT;
        else if (strcmp(arg[2],"<=") == 0) condition = LE;
        else if (strcmp(arg[2],">") == 0) condition = GT;
        else if (strcmp(arg[2],">=") == 0) condition = GE;
        else if (strcmp(arg[2],"==") == 0) condition = EQ;
        else if (strcmp(arg[2],"!=") == 0) condition = NEQ;
        else if (strcmp(arg[2],"<>") == 0) condition = BETWEEN;
        else error->all(FLERR,"Illegal group command");
        tagint bound1,bound2;
        bound1 = utils::tnumeric(FLERR,arg[3],false,lmp);
        bound2 = -1;
        if (condition == BETWEEN) {
          if (narg != 5) error->all(FLERR,"Illegal group command");
          bound2 = utils::tnumeric(FLERR,arg[4],false,lmp);
        } else if (narg != 4) error->all(FLERR,"Illegal group command");
        int *attribute = nullptr;
        tagint *tattribute = nullptr;
        if (category == TYPE) attribute = atom->type;
        else if (category == ID) tattribute = atom->tag;
        if (attribute) {
          if (condition == LT) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] < bound1) mask[i] |= bit;
          } else if (condition == LE) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] <= bound1) mask[i] |= bit;
          } else if (condition == GT) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] > bound1) mask[i] |= bit;
          } else if (condition == GE) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= bound1) mask[i] |= bit;
          } else if (condition == EQ) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] == bound1) mask[i] |= bit;
          } else if (condition == NEQ) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] != bound1) mask[i] |= bit;
          } else if (condition == BETWEEN) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= bound1 && attribute[i] <= bound2)
                mask[i] |= bit;
          }
        } else {
          if (condition == LT) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] < bound1) mask[i] |= bit;
          } else if (condition == LE) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] <= bound1) mask[i] |= bit;
          } else if (condition == GT) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] > bound1) mask[i] |= bit;
          } else if (condition == GE) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= bound1) mask[i] |= bit;
          } else if (condition == EQ) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] == bound1) mask[i] |= bit;
          } else if (condition == NEQ) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] != bound1) mask[i] |= bit;
          } else if (condition == BETWEEN) {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= bound1 && tattribute[i] <= bound2)
                mask[i] |= bit;
          }
        }
      } else {
        int *attribute = nullptr;
        tagint *tattribute = nullptr;
        if (category == TYPE) attribute = atom->type;
        else if (category == ID) tattribute = atom->tag;
        tagint start,stop,delta;
        for (int iarg = 2; iarg < narg; iarg++) {
          delta = 1;
          try {
            ValueTokenizer values(arg[iarg],":");
            start = values.next_tagint();
            if (utils::strmatch(arg[iarg],"^-?\\d+$")) {
              stop = start;
            } else if (utils::strmatch(arg[iarg],"^-?\\d+:-?\\d+$")) {
              stop = values.next_tagint();
            } else if (utils::strmatch(arg[iarg],"^-?\\d+:-?\\d+:\\d+$")) {
              stop = values.next_tagint();
              delta = values.next_tagint();
            } else throw TokenizerException("Syntax error","");
          } catch (TokenizerException &e) {
            error->all(FLERR,"Incorrect range string '{}': {}",arg[iarg],e.what());
          }
          if (delta < 1)
            error->all(FLERR,"Illegal range increment value");
          if (attribute) {
            for (i = 0; i < nlocal; i++)
              if (attribute[i] >= start && attribute[i] <= stop &&
                  (attribute[i]-start) % delta == 0) mask[i] |= bit;
          } else {
            for (i = 0; i < nlocal; i++)
              if (tattribute[i] >= start && tattribute[i] <= stop &&
                  (tattribute[i]-start) % delta == 0) mask[i] |= bit;
          }
        }
      }
    } else if (strcmp(arg[1],"variable") == 0) {
      int ivar = input->variable->find(arg[2]);
      if (ivar < 0) error->all(FLERR,"Variable name {} for group does not exist",arg[2]);
      if (!input->variable->atomstyle(ivar))
        error->all(FLERR,"Variable {} for group is invalid style",arg[2]);
      double *aflag;
      memory->create(aflag,nlocal,"group:aflag");
      input->variable->compute_atom(ivar,0,aflag,1,0);
      for (i = 0; i < nlocal; i++)
        if (aflag[i] != 0.0) mask[i] |= bit;
      memory->destroy(aflag);
    } else if (strcmp(arg[1],"include") == 0) {
      if (narg != 3) error->all(FLERR,"Illegal group include command");
      error->all(FLERR,"Unknown group include keyword {}",arg[2]);
    } else if (strcmp(arg[1],"subtract") == 0) {
      if (narg < 4) utils::missing_cmd_args(FLERR, "group subtract", error);
      int length = narg-2;
      std::vector<int> list(length);
      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR,"Group ID {} does not exist", arg[iarg]);
        if (dynamic[jgroup]) error->all(FLERR,"Cannot subtract dynamic groups");
        list[iarg-2] = jgroup;
      }
      int otherbit = bitmask[list[0]];
      for (i = 0; i < nlocal; i++)
        if (mask[i] & otherbit) mask[i] |= bit;
      int inverse = inversemask[igroup];
      for (int ilist = 1; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & otherbit) mask[i] &= inverse;
      }
    } else if (strcmp(arg[1],"union") == 0) {
      if (narg < 3) utils::missing_cmd_args(FLERR, "group union", error);
      int length = narg-2;
      std::vector<int> list(length);
      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR,"Group ID {} does not exist",arg[iarg]);
        if (dynamic[jgroup]) error->all(FLERR,"Cannot union groups from a dynamic group");
        list[iarg-2] = jgroup;
      }
      int otherbit;
      for (int ilist = 0; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & otherbit) mask[i] |= bit;
      }
    } else if (strcmp(arg[1],"intersect") == 0) {
      if (narg < 4) utils::missing_cmd_args(FLERR, "group intersect", error);
      int length = narg-2;
      std::vector<int> list(length);
      int jgroup;
      for (int iarg = 2; iarg < narg; iarg++) {
        jgroup = find(arg[iarg]);
        if (jgroup == -1) error->all(FLERR,"Group ID {} does not exist",arg[iarg]);
        if (dynamic[jgroup])
          error->all(FLERR,"Cannot intersect groups using a dynamic group");
        list[iarg-2] = jgroup;
      }
      int otherbit,ok,ilist;
      for (i = 0; i < nlocal; i++) {
        ok = 1;
        for (ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((mask[i] & otherbit) == 0) ok = 0;
        }
        if (ok) mask[i] |= bit;
      }
    } else if (strcmp(arg[1],"dynamic") == 0) {
      if (narg < 4) error->all(FLERR,"Illegal group command");
      if (strcmp(arg[0],arg[2]) == 0)
        error->all(FLERR,"Group dynamic cannot reference itself");
      if (find(arg[2]) < 0)
        error->all(FLERR,"Group dynamic parent group {} does not exist", arg[2]);
      if (igroup == 0) error->all(FLERR,"Group all cannot be made dynamic");
      if (dynamic[igroup])
        modify->delete_fix(std::string("GROUP_") + names[igroup]);
      dynamic[igroup] = 1;
      std::string fixcmd = "GROUP_";
      fixcmd += fmt::format("{} {} GROUP",names[igroup],arg[2]);
      for (i = 3; i < narg; i++) fixcmd += std::string(" ") + arg[i];
      modify->add_fix(fixcmd);
    } else if (strcmp(arg[1],"static") == 0) {
      if (narg != 2) error->all(FLERR,"Illegal group static command");
      if (dynamic[igroup])
        modify->delete_fix(std::string("GROUP_") + names[igroup]);
      dynamic[igroup] = 0;
    } else error->all(FLERR,"Unknown group command keyword: {}",arg[1]);
  } catch (LAMMPSException & e) {
    if (created) {
      delete [] names[igroup];
      names[igroup] = nullptr;
      ngroup--;
    }
    throw e;
  }
  int n;
  n = 0;
  for (i = 0; i < nlocal; i++) if (mask[i] & bit) n++;
  double rlocal = n;
  double all;
  MPI_Allreduce(&rlocal,&all,1,MPI_DOUBLE,MPI_SUM,world);
  if (me == 0) {
    if (dynamic[igroup])
      utils::logmesg(lmp,"dynamic group {} defined\n",names[igroup]);
    else
      utils::logmesg(lmp,"{:.15g} atoms in group {}\n",all,names[igroup]);
  }
}
void Group::assign(const std::string &groupcmd)
{
  auto args = utils::split_words(groupcmd);
  std::vector<char*> newarg(args.size());
  int i=0;
  for (const auto &arg : args) {
    newarg[i++] = (char *)arg.c_str();
  }
  assign(args.size(),newarg.data());
}
void Group::create(const std::string &name, int *flag)
{
  int i;
  int igroup = find(name);
  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR,"Too many groups (max {})", MAX_GROUP);
    igroup = find_unused();
    names[igroup] = utils::strdup(name);
    ngroup++;
  }
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int bit = bitmask[igroup];
  for (i = 0; i < nlocal; i++)
    if (flag[i]) mask[i] |= bit;
}
int Group::find(const std::string &name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && (name == names[igroup])) return igroup;
  return -1;
}
int Group::find_or_create(const char *name)
{
  int igroup = find(name);
  if (igroup >= 0) return igroup;
  if (ngroup == MAX_GROUP) error->all(FLERR,"Too many groups (max {})", MAX_GROUP);
  igroup = find_unused();
  names[igroup] = utils::strdup(name);
  ngroup++;
  return igroup;
}
int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == nullptr) return igroup;
  return -1;
}
void Group::write_restart(FILE *fp)
{
  fwrite(&ngroup,sizeof(int),1,fp);
  int n;
  int count = 0;
  for (int i = 0; i < MAX_GROUP; i++) {
    if (names[i]) n = strlen(names[i]) + 1;
    else n = 0;
    fwrite(&n,sizeof(int),1,fp);
    if (n) {
      fwrite(names[i],sizeof(char),n,fp);
      count++;
    }
    if (count == ngroup) break;
  }
}
void Group::read_restart(FILE *fp)
{
  int i,n;
  for (i = 0; i < MAX_GROUP; i++) delete [] names[i];
  if (me == 0) utils::sfread(FLERR,&ngroup,sizeof(int),1,fp,nullptr,error);
  MPI_Bcast(&ngroup,1,MPI_INT,0,world);
  int count = 0;
  for (i = 0; i < MAX_GROUP; i++) {
    if (count == ngroup) {
      names[i] = nullptr;
      continue;
    }
    if (me == 0) utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n) {
      names[i] = new char[n];
      if (me == 0) utils::sfread(FLERR,names[i],sizeof(char),n,fp,nullptr,error);
      MPI_Bcast(names[i],n,MPI_CHAR,0,world);
      count++;
    } else names[i] = nullptr;
  }
}
bigint Group::count_all()
{
  bigint nme = atom->nlocal;
  bigint nall;
  MPI_Allreduce(&nme,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
bigint Group::count(int igroup)
{
  int groupbit = bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;
  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
bigint Group::count(int igroup, Region *region)
{
  region->prematch();
  const int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) n++;
  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
double Group::mass(int igroup)
{
  int groupbit = bitmask[igroup];
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
double Group::mass(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
double Group::charge(int igroup)
{
  int groupbit = bitmask[igroup];
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) qone += q[i];
  double qall;
  MPI_Allreduce(&qone,&qall,1,MPI_DOUBLE,MPI_SUM,world);
  return qall;
}
double Group::charge(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
      qone += q[i];
  double qall;
  MPI_Allreduce(&qone,&qall,1,MPI_DOUBLE,MPI_SUM,world);
  return qall;
}
void Group::bounds(int igroup, double *minmax)
{
  int groupbit = bitmask[igroup];
  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }
  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];
  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);
  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}
void Group::bounds(int igroup, double *minmax, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }
  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];
  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);
  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}
void Group::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;
  double massone;
  double unwrap[3];
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }
  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::xcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;
  double massone;
  double unwrap[3];
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = rmass[i];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = mass[type[i]];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }
  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::vcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  }
  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::vcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = rmass[i];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = mass[type[i]];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  }
  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::fcm(int igroup, double *cm)
{
  int groupbit = bitmask[igroup];
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }
  MPI_Allreduce(flocal,cm,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::fcm(int igroup, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }
  MPI_Allreduce(flocal,cm,3,MPI_DOUBLE,MPI_SUM,world);
}
double Group::ke(int igroup)
{
  int groupbit = bitmask[igroup];
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}
double Group::ke(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}
double Group::gyration(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double rg = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg += (dx*dx + dy*dy + dz*dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg,&rg_all,1,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) return sqrt(rg_all/masstotal);
  return 0.0;
}
double Group::gyration(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double rg = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg += (dx*dx + dy*dy + dz*dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg,&rg_all,1,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) return sqrt(rg_all/masstotal);
  return 0.0;
}
void Group::angmom(int igroup, double *cm, double *lmom)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double p[3];
  p[0] = p[1] = p[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      p[0] += massone * (dy*v[i][2] - dz*v[i][1]);
      p[1] += massone * (dz*v[i][0] - dx*v[i][2]);
      p[2] += massone * (dx*v[i][1] - dy*v[i][0]);
    }
  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::angmom(int igroup, double *cm, double *lmom, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double p[3];
  p[0] = p[1] = p[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      p[0] += massone * (dy*v[i][2] - dz*v[i][1]);
      p[1] += massone * (dz*v[i][0] - dx*v[i][2]);
      p[2] += massone * (dx*v[i][1] - dy*v[i][0]);
    }
  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::torque(int igroup, double *cm, double *tq)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  double unwrap[3];
  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy*f[i][2] - dz*f[i][1];
      tlocal[1] += dz*f[i][0] - dx*f[i][2];
      tlocal[2] += dx*f[i][1] - dy*f[i][0];
    }
  MPI_Allreduce(tlocal,tq,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::torque(int igroup, double *cm, double *tq, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  double unwrap[3];
  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy*f[i][2] - dz*f[i][1];
      tlocal[1] += dz*f[i][0] - dx*f[i][2];
      tlocal[2] += dx*f[i][1] - dy*f[i][0];
    }
  MPI_Allreduce(tlocal,tq,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::inertia(int igroup, double *cm, double itensor[3][3])
{
  int i,j;
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];
  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}
void Group::inertia(int igroup, double *cm, double itensor[3][3], Region *region)
{
  int i,j;
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];
  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}
void Group::omega(double *angmom, double inertia[3][3], double *w)
{
  double idiag[3],ex[3],ey[3],ez[3],cross[3];
  double evectors[3][3],inverse[3][3];
  double determinant = inertia[0][0]*inertia[1][1]*inertia[2][2] +
    inertia[0][1]*inertia[1][2]*inertia[2][0] +
    inertia[0][2]*inertia[1][0]*inertia[2][1] -
    inertia[0][0]*inertia[1][2]*inertia[2][1] -
    inertia[0][1]*inertia[1][0]*inertia[2][2] -
    inertia[2][0]*inertia[1][1]*inertia[0][2];
  if (determinant > EPSILON) {
    inverse[0][0] = inertia[1][1]*inertia[2][2] - inertia[1][2]*inertia[2][1];
    inverse[0][1] = -(inertia[0][1]*inertia[2][2] -
                      inertia[0][2]*inertia[2][1]);
    inverse[0][2] = inertia[0][1]*inertia[1][2] - inertia[0][2]*inertia[1][1];
    inverse[1][0] = -(inertia[1][0]*inertia[2][2] -
                      inertia[1][2]*inertia[2][0]);
    inverse[1][1] = inertia[0][0]*inertia[2][2] - inertia[0][2]*inertia[2][0];
    inverse[1][2] = -(inertia[0][0]*inertia[1][2] -
                      inertia[0][2]*inertia[1][0]);
    inverse[2][0] = inertia[1][0]*inertia[2][1] - inertia[1][1]*inertia[2][0];
    inverse[2][1] = -(inertia[0][0]*inertia[2][1] -
                      inertia[0][1]*inertia[2][0]);
    inverse[2][2] = inertia[0][0]*inertia[1][1] - inertia[0][1]*inertia[1][0];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        inverse[i][j] /= determinant;
    w[0] = inverse[0][0]*angmom[0] + inverse[0][1]*angmom[1] +
      inverse[0][2]*angmom[2];
    w[1] = inverse[1][0]*angmom[0] + inverse[1][1]*angmom[1] +
      inverse[1][2]*angmom[2];
    w[2] = inverse[2][0]*angmom[0] + inverse[2][1]*angmom[1] +
      inverse[2][2]*angmom[2];
  } else {
    int ierror = MathEigen::jacobi3(inertia, idiag, evectors);
    if (ierror) error->all(FLERR, "Insufficient Jacobi rotations for group::omega");
    ex[0] = evectors[0][0];
    ex[1] = evectors[1][0];
    ex[2] = evectors[2][0];
    ey[0] = evectors[0][1];
    ey[1] = evectors[1][1];
    ey[2] = evectors[2][1];
    ez[0] = evectors[0][2];
    ez[1] = evectors[1][2];
    ez[2] = evectors[2][2];
    MathExtra::cross3(ex,ey,cross);
    if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);
    double max;
    max = MAX(idiag[0],idiag[1]);
    max = MAX(max,idiag[2]);
    if (idiag[0] < EPSILON*max) idiag[0] = 0.0;
    if (idiag[1] < EPSILON*max) idiag[1] = 0.0;
    if (idiag[2] < EPSILON*max) idiag[2] = 0.0;
    MathExtra::angmom_to_omega(angmom,ex,ey,ez,idiag,w);
  }
}
