#ifndef LMP_IMBALANCE_STORE_H
#define LMP_IMBALANCE_STORE_H
#include "imbalance.h"
namespace LAMMPS_NS {
class ImbalanceStore : public Imbalance {
public:
  ImbalanceStore(class LAMMPS *);
  ~ImbalanceStore() override;

public:
  int options(int, char **) override;
  void compute(double *) override;
  std::string info() override;

private:
  char *name;
};
} // namespace LAMMPS_NS
#endif
