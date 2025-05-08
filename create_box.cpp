#include <set>
#include <map>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <cstdio>
#include <string>
#include <cmath>
#include <mpi.h>
#include "lammps.h"
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "command.h"
#include "create_box.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "region_block.h"
#include "update.h"
using namespace LAMMPS_NS;
CreateBox::CreateBox(LAMMPS *lmp) : Command(lmp) {}
void CreateBox::command(int narg, char **arg) {
  domain->box_exist = 1;
  domain->triclinic = 0;
  domain->boxlo[0] = lmp->region_block->xlo;
  domain->boxhi[0] = lmp->region_block->xhi;
  domain->boxlo[1] = lmp->region_block->ylo;
  domain->boxhi[1] = lmp->region_block->yhi;
  domain->boxlo[2] = lmp->region_block->zlo;
  domain->boxhi[2] = lmp->region_block->zhi;
  atom->ntypes = utils::inumeric(FLERR, arg[0], false, lmp);
  int iarg = 2;
  update->ntimestep = 0;
  atom->allocate_type_arrays();
  atom->avec->grow(1);
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();
}
