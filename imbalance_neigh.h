#ifndef LMP_IMBALANCE_NEIGH_H
#define LMP_IMBALANCE_NEIGH_H 
#include "imbalance.h"
namespace LAMMPS_NS {
class ImbalanceNeigh : public Imbalance {
 public:
  ImbalanceNeigh(class LAMMPS *);
 public:
  int options(int, char **) override;
  void compute(double *) override;
  std::string info() override;
 private:
  double factor;
  int did_warn;
};
}
#endif
