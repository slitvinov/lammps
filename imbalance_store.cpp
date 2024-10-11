#include "imbalance_store.h"
#include "atom.h"
#include "error.h"
using namespace LAMMPS_NS;
ImbalanceStore::ImbalanceStore(LAMMPS *lmp) : Imbalance(lmp), name(nullptr) {}
ImbalanceStore::~ImbalanceStore()
{
  delete[] name;
}
int ImbalanceStore::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal balance weight command");
  name = utils::strdup(arg[0]);
  return 1;
}
void ImbalanceStore::compute(double *weight)
{
  int flag, cols;
  int index = atom->find_custom(name, flag, cols);
  if (index < 0 || flag != 1 || cols)
    error->all(FLERR, "Balance weight store vector does not exist");
  double *prop = atom->dvector[index];
  const int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; ++i) prop[i] = weight[i];
}
std::string ImbalanceStore::info()
{
  return fmt::format("  storing weight in atom property d_{}\n", name);
}
