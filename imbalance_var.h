#ifndef LMP_IMBALANCE_VAR_H
#define LMP_IMBALANCE_VAR_H
#include "imbalance.h"
namespace LAMMPS_NS {
class ImbalanceVar : public Imbalance {
public:
  ImbalanceVar(class LAMMPS *);
  ~ImbalanceVar() override;

public:
  int options(int, char **) override;
  void init(int) override;
  void compute(double *) override;
  std::string info() override;

private:
  char *name;
  int id;
};
} // namespace LAMMPS_NS
#endif
