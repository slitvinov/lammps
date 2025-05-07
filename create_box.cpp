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
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "create_box.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "region.h"
#include "update.h"
using namespace LAMMPS_NS;
CreateBox::CreateBox(LAMMPS *lmp) : Command(lmp) {}
void CreateBox::command(int narg, char **arg) {
  domain->box_exist = 1;
  auto region = domain->get_region_by_id(arg[1]);
  region->init();
  if (strcmp(region->style, "prism") != 0) {
    domain->triclinic = 0;
    domain->boxlo[0] = region->extent_xlo;
    domain->boxhi[0] = region->extent_xhi;
    domain->boxlo[1] = region->extent_ylo;
    domain->boxhi[1] = region->extent_yhi;
    domain->boxlo[2] = region->extent_zlo;
    domain->boxhi[2] = region->extent_zhi;
  }
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
