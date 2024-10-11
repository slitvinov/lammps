#include "fix_nph_sphere.h"
#include "error.h"
#include "modify.h"
using namespace LAMMPS_NS;
using namespace FixConst;
FixNPHSphere::FixNPHSphere(LAMMPS *lmp, int narg, char **arg) : FixNHSphere(lmp, narg, arg)
{
  if (tstat_flag) error->all(FLERR, "Temperature control can not be used with fix nph/sphere");
  if (!pstat_flag) error->all(FLERR, "Pressure control must be used with fix nph/sphere");
  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp/sphere", id_temp));
  tcomputeflag = 1;
  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure {}", id_press, id_temp));
  pcomputeflag = 1;
}
