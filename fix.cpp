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
int Fix::instance_total = 0;
Fix::Fix(LAMMPS *lmp, int, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr),
      vector_atom(nullptr), array_atom(nullptr),
      eatom(nullptr), vatom(nullptr), cvatom(nullptr) {
  instance_me = instance_total++;
  id = utils::strdup(arg[0]);
  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];
  style = utils::strdup(arg[2]);
  dynamic = 0;
  maxexchange = 0;
  stores_ids = 0;
  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  global_freq = -1;
  size_vector_variable = size_array_rows_variable = 0;
  comm_forward = comm_reverse = comm_border = 0;
  nevery = 1;
  global_freq = 1;
  maxeatom = maxvatom = maxcvatom = 0;
  vflag_atom = cvflag_atom = 0;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  copymode = 0;
}

