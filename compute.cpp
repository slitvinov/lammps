#include <map>
#include <unordered_set>
#include "pointers.h"
#include "compute.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include <cstring>
using namespace LAMMPS_NS;
#define DELTA 4
#define BIG MAXTAGINT
int Compute::instance_total = 0;
Compute::Compute(LAMMPS *lmp, int narg, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr), vector(nullptr),
      array(nullptr), vector_atom(nullptr), array_atom(nullptr),
      vector_local(nullptr), array_local(nullptr), extlist(nullptr),
      tlist(nullptr), vbiasall(nullptr) {
  instance_me = instance_total++;
  if (narg < 3)
    error->all(FLERR, "Illegal compute command");
  id = utils::strdup(arg[0]);
  if (!utils::is_id(id))
    error->all(FLERR,
               "Compute ID must be alphanumeric or underscore characters");
  igroup = group->find(arg[1]);
  if (igroup == -1)
    error->all(FLERR, "Could not find compute group ID");
  groupbit = group->bitmask[igroup];
  style = utils::strdup(arg[2]);
  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = pergrid_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;
  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  create_attribute = 0;
  tempbias = 0;
  timeflag = 0;
  comm_forward = comm_reverse = 0;
  dynamic = 0;
  dynamic_group_allow = 1;
  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_peratom = invoked_local = -1;
  invoked_flag = INVOKED_NONE;
  extra_dof = domain->dimension;
  dynamic_user = 0;
  fix_dof = 0;
  ntime = maxtime = 0;
  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  copymode = 0;
}
Compute::~Compute() {
  if (copymode)
    return;
  delete[] id;
  delete[] style;
  memory->destroy(tlist);
}
void Compute::modify_params(int narg, char **arg) {
  if (narg == 0)
    error->all(FLERR, "Illegal compute_modify command");
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "extra") == 0 ||
        strcmp(arg[iarg], "extra/dof") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute_modify command");
      extra_dof = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dynamic") == 0 ||
               strcmp(arg[iarg], "dynamic/dof") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute_modify command");
      dynamic_user = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute_modify command");
  }
}
void Compute::adjust_dof_fix() {
  fix_dof = 0;
  for (auto &ifix : modify->get_fix_list())
    if (ifix->dof_flag)
      fix_dof += ifix->dof(igroup);
}
void Compute::reset_extra_dof() { extra_dof = domain->dimension; }
void Compute::reset_extra_compute_fix(const char *) {
  error->all(FLERR,
             "Compute does not allow an extra compute or fix to be reset");
}
void Compute::addstep(bigint ntimestep) {
  int i;
  for (i = ntime - 1; i >= 0; i--) {
    if (ntimestep == tlist[i])
      return;
    if (ntimestep < tlist[i])
      break;
  }
  i++;
  if (ntime == maxtime) {
    maxtime += DELTA;
    memory->grow(tlist, maxtime, "compute:tlist");
  }
  for (int j = ntime - 1; j >= i; j--)
    tlist[j + 1] = tlist[j];
  tlist[i] = ntimestep;
  ntime++;
}
int Compute::matchstep(bigint ntimestep) {
  for (int i = ntime - 1; i >= 0; i--) {
    if (ntimestep < tlist[i])
      return 0;
    if (ntimestep == tlist[i])
      return 1;
    if (ntimestep > tlist[i])
      ntime--;
  }
  return 0;
}
void Compute::clearstep() { ntime = 0; }
