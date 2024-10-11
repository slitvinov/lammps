#include "update.h"
#include "style_integrate.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "neighbor.h"
#include <cstring>
using namespace LAMMPS_NS;
template <typename T> static Integrate *integrate_creator(LAMMPS *lmp, int narg, char **arg)
{
  return new T(lmp, narg, arg);
}
Update::Update(LAMMPS *lmp) : Pointers(lmp)
{
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
#define INTEGRATE_CLASS 
#define IntegrateStyle(key,Class) (*integrate_map)[#key] = &integrate_creator<Class>;
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS
  str = (char *) "verlet";
  create_integrate(1, &str, 1);
}
Update::~Update()
{
  delete[] unit_style;
  delete[] integrate_style;
  delete integrate;
  delete integrate_map;
}
void Update::init()
{
  if (whichflag == 0) return;
  if (whichflag == 1)
    integrate->init();
  first_update = 1;
}
void Update::set_units(const char *style)
{
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
  } else if (strcmp(style, "real") == 0) {
    force->boltz = 0.0019872067;
    force->hplanck = 95.306976368;
    force->mvv2e = 48.88821291 * 48.88821291;
    force->ftm2v = 1.0 / 48.88821291 / 48.88821291;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;
    force->qe2f = 23.060549;
    force->vxmu2f = 1.4393264316e4;
    force->xxt2kmu = 0.1;
    force->e_mass = 1.0 / 1836.1527556560675;
    force->hhmrr2e = 0.0957018663603261;
    force->mvh2r = 1.5339009481951;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;
    dt = 1.0;
    neighbor->skin = 2.0;
  } else if (strcmp(style, "metal") == 0) {
    force->boltz = 8.617343e-5;
    force->hplanck = 4.135667403e-3;
    force->mvv2e = 1.0364269e-4;
    force->ftm2v = 1.0 / 1.0364269e-4;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 1.6021765e6;
    force->qqr2e = 14.399645;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0e-3;
    force->qelectron = 1.0;
    dt = 0.001;
    neighbor->skin = 2.0;
  } else if (strcmp(style, "si") == 0) {
    force->boltz = 1.3806504e-23;
    force->hplanck = 6.62606896e-34;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e9;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-10;
    force->femtosecond = 1.0e-15;
    force->qelectron = 1.6021765e-19;
    dt = 1.0e-8;
    neighbor->skin = 0.001;
  } else if (strcmp(style, "cgs") == 0) {
    force->boltz = 1.3806504e-16;
    force->hplanck = 6.62606896e-27;
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
    force->angstrom = 1.0e-8;
    force->femtosecond = 1.0e-15;
    force->qelectron = 4.8032044e-10;
    dt = 1.0e-8;
    neighbor->skin = 0.1;
  } else if (strcmp(style, "electron") == 0) {
    force->boltz = 3.16681534e-6;
    force->hplanck = 0.1519829846;
    force->mvv2e = 1.06657236;
    force->ftm2v = 0.937582899;
    force->mv2d = 1.0;
    force->nktv2p = 2.94210108e13;
    force->qqr2e = 1.0;
    force->qe2f = 1.94469051e-10;
    force->vxmu2f = 3.39893149e1;
    force->xxt2kmu = 3.13796367e-2;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.88972612;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;
    dt = 0.001;
    neighbor->skin = 2.0;
  } else if (strcmp(style, "micro") == 0) {
    force->boltz = 1.3806504e-8;
    force->hplanck = 6.62606896e-13;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.987556e6;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-4;
    force->femtosecond = 1.0e-9;
    force->qelectron = 1.6021765e-7;
    dt = 2.0;
    neighbor->skin = 0.1;
  } else if (strcmp(style, "nano") == 0) {
    force->boltz = 0.013806504;
    force->hplanck = 6.62606896e-4;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 230.7078669;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-1;
    force->femtosecond = 1.0e-6;
    force->qelectron = 1.0;
    dt = 0.00045;
    neighbor->skin = 0.1;
  } else
    error->all(FLERR, "Illegal units command");
  delete[] unit_style;
  unit_style = utils::strdup(style);
  if (!dt_default && (comm->me == 0)) {
    error->warning(FLERR, "Changing timestep from {:.6} to {:.6} due to changing units to {}",
                   dt_old, dt, unit_style);
  }
  dt_default = 1;
}
void Update::create_integrate(int narg, char **arg, int trysuffix)
{
  if (narg < 1) error->all(FLERR, "Illegal run_style command");
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
void Update::new_integrate(char *style, int narg, char **arg)
{
  if (integrate_map->find(style) != integrate_map->end()) {
    IntegrateCreator &integrate_creator = (*integrate_map)[style];
    integrate = integrate_creator(lmp, narg, arg);
    return;
  }
  error->all(FLERR, "Illegal integrate style");
}
void Update::reset_timestep(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "reset_timestep", error);
  reset_timestep(utils::bnumeric(FLERR, arg[0], false, lmp), true);
  if (narg > 1) {
    if (strcmp(arg[1], "time") == 0) {
      if (narg < 3) utils::missing_cmd_args(FLERR, "reset_timestep time", error);
      atimestep = ntimestep;
      atime = utils::numeric(FLERR, arg[2], false, lmp);
    } else
      error->all(FLERR, "Unknown reset_timestep option {}", arg[1]);
  }
}
void Update::reset_timestep(bigint newstep, bool do_check)
{
  if (newstep < 0) error->all(FLERR, "Timestep must be >= 0");
  bigint oldstep = ntimestep;
  ntimestep = newstep;
  if (newstep >= oldstep) update_time();
  if (newstep < oldstep) {
    atime = 0.0;
    atimestep = newstep;
  }
  if (do_check) {
    for (const auto &ifix : modify->get_fix_list())
      if (ifix->time_depend)
        error->all(FLERR, "Cannot reset timestep with time-dependent fix {} defined", ifix->style);
  }
  eflag_global = vflag_global = -1;
  for (const auto &icompute : modify->get_compute_list()) {
    icompute->invoked_scalar = -1;
    icompute->invoked_vector = -1;
    icompute->invoked_array = -1;
    icompute->invoked_peratom = -1;
    icompute->invoked_local = -1;
    if (icompute->timeflag) icompute->clearstep();
  }
  neighbor->reset_timestep(ntimestep);
}
void Update::update_time()
{
  atime += (ntimestep - atimestep) * dt;
  atimestep = ntimestep;
}
double Update::memory_usage()
{
  double bytes = 0;
  if (whichflag == 1)
    bytes += integrate->memory_usage();
  return bytes;
}
