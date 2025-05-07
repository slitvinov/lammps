#include <vector>
#include <cstring>
#include <map>
#include <cstdio>
#include <string>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "utils.h"
#include "update.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "neighbor.h"
#include "lmptype.h"
#include "integrate.h"
#include "verlet.h"
using namespace LAMMPS_NS;
template <typename T>
static Integrate *integrate_creator(LAMMPS *lmp, int narg, char **arg) {
  return new T(lmp, narg, arg);
}
Update::Update(LAMMPS *lmp) : Pointers(lmp) {
  char *str;
  ntimestep = 0;
  atime = 0.0;
  atimestep = 0;
  first_update = 0;
  whichflag = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  restrict_output = 0;
  setupflag = 0;
  multireplica = 0;
  eflag_global = vflag_global = -1;
  eflag_atom = vflag_atom = 0;
  dt_default = 1;
  dt = 0.0;
  unit_style = nullptr;
  set_units("lj");
  integrate_style = nullptr;
  integrate = nullptr;
  integrate_map = new IntegrateCreatorMap();
  (*integrate_map)["verlet"] = &integrate_creator<Verlet>;
  str = (char *)"verlet";
  create_integrate(1, &str, 1);
}
void Update::init() {
  if (whichflag == 0)
    return;
  if (whichflag == 1)
    integrate->init();
  first_update = 1;
}
void Update::set_units(const char *style) {
  double dt_old = dt;
  if (strcmp(style, "lj") == 0) {
    force->boltz = 1.0;
    force->hplanck = 1.0;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;
    dt = 0.005;
    neighbor->skin = 0.3;
  }
  delete[] unit_style;
  unit_style = utils::strdup(style);
  dt_default = 1;
}
void Update::create_integrate(int narg, char **arg, int trysuffix) {
  delete[] integrate_style;
  delete integrate;
  if (narg - 1 > 0) {
    new_integrate(arg[0], narg - 1, &arg[1]);
  } else {
    new_integrate(arg[0], 0, nullptr);
  }
  std::string estyle = arg[0];
  integrate_style = utils::strdup(estyle);
}
void Update::new_integrate(char *style, int narg, char **arg) {
  if (integrate_map->find(style) != integrate_map->end()) {
    IntegrateCreator &integrate_creator = (*integrate_map)[style];
    integrate = integrate_creator(lmp, narg, arg);
    return;
  }
}
void Update::update_time() {
  atime += (ntimestep - atimestep) * dt;
  atimestep = ntimestep;
}
