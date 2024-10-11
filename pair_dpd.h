#ifdef PAIR_CLASS
PairStyle(dpd,PairDPD);
#else
#ifndef LMP_PAIR_DPD_H
#define LMP_PAIR_DPD_H 
#include "pair.h"
namespace LAMMPS_NS {
class PairDPD : public Pair {
 public:
  PairDPD(class LAMMPS *);
  ~PairDPD() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
 protected:
  double cut_global, temperature;
  double special_sqrt[4];
  int seed;
  double **cut;
  double **a0, **gamma;
  double **sigma;
  class RanMars *random;
  virtual void allocate();
};
}
#endif
#endif
