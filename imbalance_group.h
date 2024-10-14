#ifndef LMP_IMBALANCE_GROUP_H
#define LMP_IMBALANCE_GROUP_H
#include "imbalance.h"
namespace LAMMPS_NS {
class ImbalanceGroup : public Imbalance {
public:
  ImbalanceGroup(class LAMMPS *);
  ~ImbalanceGroup() override;
  int options(int, char **) override;
  void compute(double *) override;
  std::string info() override;

private:
  int num;
  int *id;
  double *factor;
};
} // namespace LAMMPS_NS
#endif
