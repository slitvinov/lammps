#include <set>
#include <map>
#include <vector>
#include <string>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "fix.h"
#include "atom.h"
#include "atom_masks.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include <cstring>
using namespace LAMMPS_NS;
using namespace FixConst;
int Fix::instance_total = 0;
Fix::Fix(LAMMPS *lmp, int, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr), extlist(nullptr),
      vector_atom(nullptr), array_atom(nullptr), vector_local(nullptr),
      array_local(nullptr), eatom(nullptr), vatom(nullptr), cvatom(nullptr) {
  instance_me = instance_total++;
  id = utils::strdup(arg[0]);
  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];
  style = utils::strdup(arg[2]);
  restart_global = restart_peratom = restart_file = 0;
  force_reneighbor = 0;
  box_change = NO_BOX_CHANGE;
  thermo_energy = 0;
  thermo_virial = 0;
  energy_global_flag = energy_peratom_flag = 0;
  virial_global_flag = virial_peratom_flag = 0;
  ecouple_flag = 0;
  rigid_flag = 0;
  no_change_box = 0;
  time_integrate = 0;
  time_depend = 0;
  create_attribute = 0;
  restart_pbc = 0;
  wd_header = wd_section = 0;
  dynamic = 0;
  dof_flag = 0;
  special_alter_flag = 0;
  enforce2d_flag = 0;
  respa_level_support = 0;
  respa_level = -1;
  maxexchange = 0;
  maxexchange_dynamic = 0;
  pre_exchange_migrate = 0;
  stores_ids = 0;
  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = pergrid_flag = 0;
  global_freq = local_freq = peratom_freq = pergrid_freq = -1;
  size_vector_variable = size_array_rows_variable = 0;
  comm_forward = comm_reverse = comm_border = 0;
  restart_reset = 0;
  nevery = 1;
  global_freq = 1;
  maxeatom = maxvatom = maxcvatom = 0;
  vflag_atom = cvflag_atom = 0;
  centroidstressflag = CENTROID_SAME;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  forward_comm_device = 0;
  copymode = 0;
}

