#include "imbalance_group.h"
#include "atom.h"
#include "error.h"
#include "group.h"
using namespace LAMMPS_NS;
ImbalanceGroup::ImbalanceGroup(LAMMPS *lmp) : Imbalance(lmp), id(nullptr), factor(nullptr) {}
ImbalanceGroup::~ImbalanceGroup()
{
  delete[] id;
  delete[] factor;
}
int ImbalanceGroup::options(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR, "Illegal balance weight command");
  num = utils::inumeric(FLERR, arg[0], false, lmp);
  if (num < 1) error->all(FLERR, "Illegal balance weight command");
  if (2 * num + 1 > narg) error->all(FLERR, "Illegal balance weight command");
  id = new int[num];
  factor = new double[num];
  for (int i = 0; i < num; ++i) {
    id[i] = group->find(arg[2 * i + 1]);
    if (id[i] < 0) error->all(FLERR, "Unknown group in balance weight command: {}", arg[2 * i + 1]);
    factor[i] = utils::numeric(FLERR, arg[2 * i + 2], false, lmp);
    if (factor[i] <= 0.0) error->all(FLERR, "Illegal balance weight command");
  }
  return 2 * num + 1;
}
void ImbalanceGroup::compute(double *weight)
{
  const int *const mask = atom->mask;
  const int *const bitmask = group->bitmask;
  const int nlocal = atom->nlocal;
  if (num == 0) return;
  for (int i = 0; i < nlocal; ++i) {
    const int imask = mask[i];
    for (int j = 0; j < num; ++j) {
      if (imask & bitmask[id[j]]) weight[i] *= factor[j];
    }
  }
}
std::string ImbalanceGroup::info()
{
  std::string mesg;
  if (num > 0) {
    const char *const *const names = group->names;
    mesg += "  group weights:";
    for (int i = 0; i < num; ++i) mesg += fmt::format(" {}={}", names[id[i]], factor[i]);
    mesg += "\n";
  }
  return mesg;
}
