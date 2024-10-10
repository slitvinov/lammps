/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "create_box.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "region.h"
#include "region_prism.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "create_box", error);

  if (domain->box_exist) error->all(FLERR, "Cannot create_box after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR, "Cannot run 2d simulation with nonperiodic Z dimension");

  domain->box_exist = 1;

  // region check

  auto region = domain->get_region_by_id(arg[1]);
  if (!region) error->all(FLERR, "Create_box region {} does not exist", arg[1]);
  if (region->bboxflag == 0) error->all(FLERR, "Create_box region does not support a bounding box");

  region->init();

  // if region not prism:
  //   setup orthogonal domain
  //   set simulation domain from region extent
  // if region is prism:
  //   seutp triclinic domain
  //   set simulation domain params from prism params

  if (strcmp(region->style, "prism") != 0) {
    domain->triclinic = 0;
    domain->boxlo[0] = region->extent_xlo;
    domain->boxhi[0] = region->extent_xhi;
    domain->boxlo[1] = region->extent_ylo;
    domain->boxhi[1] = region->extent_yhi;
    domain->boxlo[2] = region->extent_zlo;
    domain->boxhi[2] = region->extent_zhi;

  } else {
    domain->triclinic = 1;
    auto prism = dynamic_cast<RegPrism *>(region);
    domain->boxlo[0] = prism->xlo;
    domain->boxhi[0] = prism->xhi;
    domain->boxlo[1] = prism->ylo;
    domain->boxhi[1] = prism->yhi;
    domain->boxlo[2] = prism->zlo;
    domain->boxhi[2] = prism->zhi;
    domain->xy = prism->xy;
    domain->xz = prism->xz;
    domain->yz = prism->yz;
  }

  // set atom and topology type quantities

  atom->ntypes = utils::inumeric(FLERR, arg[0], false, lmp);

  // process optional args that can overwrite default settings
  int iarg = 2;

  // problem setup using info from header
  // deallocate/grow ensures any extra settings are used for topology arrays
  // necessary in case no create_atoms is performed

  update->ntimestep = 0;

  atom->allocate_type_arrays();
  atom->avec->grow(1);

  domain->print_box("Created ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();
}
